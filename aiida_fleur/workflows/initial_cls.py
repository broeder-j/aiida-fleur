# -*- coding: utf-8 -*-
###############################################################################
# Copyright (c), Forschungszentrum Jülich GmbH, IAS-1/PGI-1, Germany.         #
#                All rights reserved.                                         #
# This file is part of the AiiDA-FLEUR package.                               #
#                                                                             #
# The code is hosted on GitHub at https://github.com/broeder-j/aiida-fleur    #
# For further information on the license, see the LICENSE.txt file            #
# For further information please visit http://www.flapw.de or                 #
# http://aiida-fleur.readthedocs.io/en/develop/                               #
###############################################################################

"""
This is the worklfow 'initial_cls' using the Fleur code calculating
corelevel shifts with different methods.
"""
#TODO parsing of eigenvalues of LOS!
#TODO error handling of scf
#TODO Check if calculations failed, and termine the workflow without a raised execption
# currently the result extraction part will fail if calculations failed
#TODO USE SAME PARAMETERS! (maybe extract method for fleurinp needed)
# TODO: Allow for providing referenes as scf_ouputparameter nodes
# TODO: maybe launch all scfs at the same time
# TODO: gives only a warning currently if ref not found.
# but should lead to error if no ref is found for what should be calculated
from string import digits
from aiida.work.run import submit
from aiida.work.workchain import ToContext, WorkChain, if_
from aiida.work import workfunction as wf
from aiida.orm import Code, DataFactory, CalculationFactory, load_node, Group
from aiida.orm.querybuilder import QueryBuilder
from aiida.common.exceptions import NotExistent
from aiida_fleur.calculation.fleur import FleurCalculation
from aiida_fleur.workflows.scf import fleur_scf_wc
from aiida_fleur.tools.common_fleur_wf_util import get_natoms_element


StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
RemoteData = DataFactory('remote')
FleurinpData = DataFactory('fleur.fleurinp')
FleurProcess = FleurCalculation.process()
FleurCalc = CalculationFactory('fleur.fleur')


class fleur_initial_cls_wc(WorkChain):
    """
    Turn key solution for the calculation of core level shift


    'method' : ['initial', 'full_valence ch', 'half_valence_ch', 'ch', ...]
    'Bes' : [W4f, Be1s]
    'CLS' : [W4f, Be1s]
 toms' : ['all', 'postions' : []]
    #'references' : ['calculate', and use # calculate : 'all' , or 'calculate' : ['W', 'Be']
    'references' : { 'W': [calc/ouputnode or  fleurinp, or structure data or
                     structure data + Parameter  ], 'Be' : }
    'scf_para' : {...}, 'default'
    'relax' : True
    'relax_mode': ['Fleur', 'QE Fleur', 'QE']
    'relax_para' : {...}, 'default'
    'calculate_doses' : False
    'dos_para' : {...}, 'default'

    # defaults
    default wf_Parameters::
    'method' : 'initial'
    'atoms' : 'all
    'references' : 'calculate'
    'scf_para' : 'default'
    'relax' : True
    'relax_mode': 'QE Fleur'
    'relax_para' : 'default'
    'calculate_doses' : False
    'dos_para' : 'default'
    """

    _workflowversion = "0.3.4"
    _default_wf_para = {'references' : {},
                       'relax' : True,
                       'relax_mode': 'Fleur',
                       'relax_para' : 'default',
                       'scf_para' : 'default',
                       'same_para' : True,
                       'serial' : False}

    _default_options = {'resources' : {"num_machines": 1},
                        'max_wallclock_seconds' : 6*60*60,
                        'queue_name' : '',
                        'custom_scheduler_commands' : '',
                        #'max_memory_kb' : None,
                        'import_sys_environment' : False,
                        'environment_variables' : {}}
    
                        
    ERROR_INVALID_INPUT_RESOURCES = 1
    ERROR_INVALID_INPUT_RESOURCES_UNDERSPECIFIED = 2
    ERROR_INVALID_CODE_PROVIDED = 3
    ERROR_INPGEN_CALCULATION_FAILED = 4
    ERROR_CHANGING_FLEURINPUT_FAILED = 5
    ERROR_CALCULATION_INVALID_INPUT_FILE = 6
    ERROR_FLEUR_CALCULATION_FALIED = 7
    ERROR_CONVERGENCE_NOT_ARCHIVED = 8
    ERROR_REFERENCE_MISSING = 9
    

    @classmethod
    def define(cls, spec):
        super(fleur_initial_cls_wc, cls).define(spec)
        spec.input("wf_parameters", valid_type=ParameterData, required=False,
                   default=ParameterData(dict=cls._default_wf_para))
        spec.input("fleurinp", valid_type=FleurinpData, required=False)
        spec.input("fleur", valid_type=Code, required=True)
        spec.input("inpgen", valid_type=Code, required=False)
        spec.input("structure", valid_type=StructureData, required=False)
        spec.input("calc_parameters", valid_type=ParameterData, required=False)
        spec.input("options", valid_type=ParameterData, required=False, 
                   default=ParameterData(dict=cls._default_options))

        spec.outline(
            cls.check_input,
            cls.get_references,
            cls.run_fleur_scfs,
            if_(cls.relaxation_needed)(
                cls.relax),
            cls.find_parameters,
            cls.run_scfs_ref,
            cls.return_results
        )


    def check_input(self):
        """
        Init same context and check what input is given if it makes sence
        """
        ### input check ### ? or done automaticly, how optional?

        msg = ("INFO: Started initial_state_CLS workflow version {} "
               "Workchain node identifiers: "#{}"
               "".format(self._workflowversion))#, ProcessRegistry().current_calc_node))
        self.report(msg)

        # init
        self.ctx.last_calc = None
        self.ctx.eximated_jobs = 0
        self.ctx.run_jobs = 0
        self.ctx.calcs_res = []
        self.ctx.labels = []
        self.ctx.ref_labels = []
        self.ctx.calcs_torun = []
        self.ctx.ref_calcs_torun = []
        self.ctx.ref_calcs_res = []
        self.ctx.struc_to_relax = []
        self.ctx.successful = True
        self.ctx.warnings = []
        self.ctx.errors = []
        self.ctx.ref = {}
        self.ctx.calculate_formation_energy = True

        #Style: {atomtype : listof all corelevel, atomtype_coresetup... }
        #ie: { 'W-1' : [shift_1s, ... shift 7/2 4f],
        #      'W-1_coreconfig' : ['1s','2s',...],
        #      'W-2' : [...], 'Be-1': [], ...} #all in eV!
        self.ctx.CLS = {}
        self.ctx.cl_energies = {}# same style as CLS only energy <-> shift
        self.ctx.ref_cl_energies = {}
        #Style: {'Compound' : energy, 'ref_x' : energy , ...}
        #i.e {'Be12W' : 0.0, 'Be' : 0.104*htr_eV , 'W' : 0.12*htr_eV} # all in eV!
        self.ctx.fermi_energies = {}
        self.ctx.bandgaps = {}
        self.ctx.atomtypes = {}
        # set values, or defaults for Wf_para
        wf_dict = self.inputs.wf_parameters.get_dict()
        default = self._default_wf_para

        self.ctx.serial = wf_dict.get('serial', default.get('serial'))
        self.ctx.same_para = wf_dict.get('same_para', default.get('same_para'))
        self.ctx.scf_para = wf_dict.get('scf_para', default.get('scf_para'))
        self.ctx.relax = wf_dict.get('relax', default.get('relax'))
        self.ctx.relax_mode = wf_dict.get('relax_mode', default.get('relax_mode'))
        self.ctx.relax_para = wf_dict.get('relax_para', default.get('dos_para'))
        self.ctx.ref_calcs_res_given = []
        
        defaultoptions = self._default_options
        if self.inputs.options:
            options = self.inputs.options.get_dict()
        else:
            options = defaultoptions
        for key, val in defaultoptions.iteritems():
            options[key] = options.get(key, val)
        self.ctx.options = options

        # check if inputs given make sense # TODO sort this out in common wc
        inputs = self.inputs
        if 'fleurinp' in inputs:
            #TODO make a check if an extracted structure exists, since get_structuredata is wf
            structure = inputs.fleurinp.get_structuredata(inputs.fleurinp)
            self.ctx.elements = list(structure.get_composition().keys())
            self.ctx.calcs_torun.append(inputs.get('fleurinp'))
            #print('here1')
            if 'structure' in inputs:
                warning = 'WARNING: Ignoring Structure input, because Fleurinp was given'
                self.ctx.warnings.append(warning)
                self.report(warning)
            if 'calc_parameters' in inputs:
                warning = 'WARNING: Ignoring parameter input, because Fleurinp was given'
                self.ctx.warnings.append(warning)
                self.report(warning)
        elif 'structure' in inputs:
            self.ctx.elements = list(inputs.structure.get_composition().keys())
            #self.ctx.elements = list(s.get_symbols_set())
            if 'inpgen' not in inputs:
                error = 'ERROR: StructureData was provided, but no inpgen code was provided'
                self.ctx.errors.append(error)
                self.control_end_wc(error)
                return self.ERROR_INVALID_INPUT_RESOURCES
            if 'calc_parameters' in inputs:
                self.ctx.calcs_torun.append(
                    [inputs.get('structure'), inputs.get('calc_parameters')])
            else:
                self.ctx.calcs_torun.append(inputs.get('structure'))
        else:
            error = 'ERROR: No StructureData nor FleurinpData was provided'
            self.ctx.errors.append(error)
            self.control_end_wc(error)
            return self.ERROR_INVALID_INPUT_RESOURCES
        self.report('INFO: elements in structure: {}'.format(self.ctx.elements))


    def get_references(self):
        """
        To calculate a CLS in initial state approx, we need reference calculations
        to the Elemental crystals. First it is checked if the user has provided them
        Second the database is checked, if there are structures with certain extras.
        Third the COD database is searched for the elemental Cystal structures.
        If some referneces are not found stop here.
        Are there already calculation of these 'references', ggf use them.
        We do not put these calculation in the calculation queue yet because we
        need specific parameters for them
        """

        self.report('INFO: In Get_references initial_state_CLS workflow')

        references = self.inputs.wf_parameters.get_dict().get('references', {})
        # should be of the form of
        #'references' : { 'W': calc, outputnode of workflow or fleurinp,
        #                 or structure data or (structure data + Parameter),
        #                 'Be' : ...}
        self.ctx.ref_calcs_torun = []
        self.ctx.ref = {}
        self.ctx.abort = False

        struc_group = references.get('group', None)
        para_group = references.get('para_group', None)

        #TODO better checks if ref makes sense?

        # get specific element reference if given override
        #print(self.ctx.elements)
        elements = self.ctx.elements # ggf copy because ctx.elements will be modified
        for elem in elements:
            #to_calc[elem] = 'find'
            ref_el = references.get(elem, None)
            #print ref_el
            if ref_el:
                # loading nodes
                if isinstance(ref_el, list):
                    ref_el_node = []
                    for ref_el_el in ref_el:
                        try:
                            ref_el_nodes = load_node(ref_el_el)
                        except:
                            ref_el_node = None
                            self.report('ERROR: The reference node in the list '
                                        '(id or uuid) provided: {} for element: '
                                        '{} could not be loaded with load_node'
                                        ''.format(ref_el_el, elem))
                            self.ctx.abort = True
                        ref_el_node.append(ref_el_nodes)
                else:
                    try:
                        ref_el_node = load_node(ref_el)
                    except:# NotExistent: No node was found
                        ref_el_node = None
                        self.report('ERROR: The reference node (id or uuid) '
                                    'provided: {} for element: {} could'
                                    'not be loaded with load_node'
                                    ''.format(ref_el, elem))
                        self.ctx.abort = True

                # expecting nodes and filling ref_calcs_torun
                if isinstance(ref_el_node, list):#(StructureData, ParameterData)):
                    #enforced parameters, add directly to run queue
                    # TODO: if a scf with these parameters was already done link to it
                    # and extract the results instead of running the calculation again....
                    if len(ref_el_node) == 2:
                        if isinstance(ref_el_node[0], StructureData) and isinstance(ref_el_node[1], ParameterData):
                            self.ctx.ref_calcs_torun.append(ref_el_node)
                        else:
                            self.report('WARNING: I did not undestand the list with length 2 '
                                  'you gave me as reference input')
                    else:
                        self.report('WARNING: I did not undestand the list {} with length {} '
                              'you gave me as reference input'
                              ''.format(ref_el_node, len(ref_el_node)))
                elif isinstance(ref_el_node, FleurCalc):
                    #extract from fleur calc TODO
                    self.ctx.ref_cl_energies[elem] = {}
                elif isinstance(ref_el_node, ParameterData):
                    #extract from workflow output TODO
                    self.ctx.ref_cl_energies[elem] = {}
                    # test if scf_output_para, or fleur_outputpara
                    #self.ctx.ref_calcs_res_given = ref_given
                    #ref_given = []
                    #uuid_fleur_calc = 
                    #ref_given.append(uuid_fleur_calc)
                elif isinstance(ref_el_node, FleurinpData):
                    # add to calculations
                    #enforced parameters, add directly to run queue
                    self.ctx.ref_calcs_torun.append(ref_el_node)
                    #self.ctx.ref[elem] = ref_el
                elif isinstance(ref_el_node, StructureData):
                    self.ctx.ref[elem] = ref_el_node
                    self.ctx.ref_calcs_torun.append(ref_el_node)
                #elif isinstance(ref_el, initial_state_CLS):
                #    extract TODO
                else:
                    error = ("ERROR: I do not know what to do with this given "
                             "reference {} for element {}".format(ref_el, elem))
                    #print(error)
                    self.report(error)
                    self.ctx.errors.append(error)
                    self.ctx.abort = True
            elif struc_group:
                #print('here, looking in group')
                #print(elem, struc_group)
                structure, report = get_ref_from_group(elem, struc_group)
                if report:
                    self.report(report)
                parameter, report = get_para_from_group(elem, para_group)
                if structure and parameter:
                    self.ctx.ref[elem] = structure
                    self.ctx.ref_calcs_torun.append([structure, parameter])
                elif structure:
                    self.ctx.ref[elem] = structure
                    self.ctx.ref_calcs_torun.append(structure)
                else:
                    pass # report not found?


            #elif query_for_ref: # no ref given, we have to look for it.
            #    structure = querry_for_ref_structure(elem)
            #    #TODO: Advance this querry, if a calculation with the given
            #    #parameter was already done use these results
            #    if structure:
            #        self.ctx.ref[elem] = structure
            #        self.ctx.ref_calcs_torun.append(structure)# tempoary later check parameters
            #    else: #not found
            #        error = ("ERROR: Reference structure for element: {} not found."
            #                 "checkout the 'querry_for_ref_structure' method."
            #                 "to see what extras are querried for.".format(elem))
            #        #print(error)
            #        self.ctx.errors.append(error)
            #        self.ctx.abort = True
            #        self.report(error)

            else: # no reference for element found
                # do we not want to calculate it or is this an error?
                warning = ('WARNING: I did not find a reference for element {}. '
                           'If you do not calculate shifts for this element '
                           'ignore this warning. If I should calculate this will '
                           'lead to an error later. Note: without all references'
                           ' I cannot calculte the binding energy'.format(elem))
                self.report(warning)
                self.ctx.calculate_formation_energy = False
                # delete element from element list (no calculations will be launched for it)
                valid_elm = self.ctx.elements
                i = valid_elm.index(elem)
                del valid_elm[i]
                self.ctx.elements = valid_elm

        if self.ctx.abort:
            error = ('ERROR: Something was wrong with the reference input provided. '
                     'I cannot calculate from the input, or what I have found '
                     'what you want me to do. Please check the workchain report'
                     'for details.')
            self.control_end_wc(error)
            return self.ERROR_REFERENCE_MISSING

        self.report('INFO ref_calcs_torun: {} '.format(self.ctx.ref_calcs_torun))

        # check if a structureData for these elements was given
        #if yes add to ref_calc to run
        #was also a prameter node given for the element?
        #yes run with these
        #no was on given for the host structure, extract element parameternode

        #else use parameters extracted from host calculation # TODO

        #check if there is a structure from this element in the database with extras:
        # with extra.type = 'bulk', extra.specific = 'reference',
        #'extra.elemental' = True, extra.structure = 'W'
        # check if input parameter node values for the calculation are the same.

        #if yes, if a calculation exists use that result
        #else do a calculation on that structure as above



    def run_fleur_scfs(self):
        """
        Run SCF-cycles for all structures, calculations given in certain workflow arrays.
        """
        self.report('INFO: In run_fleur_scfs initial_state_CLS workflow')

        para = self.ctx.scf_para
        if para == 'default':
            wf_parameter = {}
        else:
            wf_parameter = para
        wf_parameter['serial'] = self.ctx.serial
        #wf_parameter['options'] = self.ctx.options        
        wf_parameters = ParameterData(dict=wf_parameter)
        res_all = []
        options = ParameterData(dict=self.ctx.options)
        # for each calulation in self.ctx.calcs_torun #TODO what about wf params?
        res = None
        #print(self.ctx.calcs_torun)
        for node in self.ctx.calcs_torun:
            #print node
            scf_label = 'cls|scf_wc main'
            scf_description = 'cls|scf of the main structure'
            if isinstance(node, StructureData):
                res = self.submit(fleur_scf_wc, wf_parameters=wf_parameters, structure=node,
                             inpgen=self.inputs.inpgen, fleur=self.inputs.fleur, options=options,
                             label=scf_label, description=scf_description)#
            #elif isinstance(node, FleurinpData):
            #    res = fleur_scf_wc.run(wf_parameters=wf_parameters, structure=node,
            #                inpgen = self.inputs.inpgen, fleur=self.inputs.fleur)#
            elif isinstance(node, list):#(StructureData, ParameterData)):
                if len(node) == 2:
                    res = self.submit(fleur_scf_wc, wf_parameters=wf_parameters,
                                 structure=node[0], calc_parameters=node[1], options=options,
                                 inpgen=self.inputs.inpgen, fleur=self.inputs.fleur,
                                 label=scf_label, description=scf_description)
                else:
                    self.report('ERROR: something in calcs_torun which I do not'
                                'recognize, list has not 2 entries: {}'.format(node))
            else:
                self.report('ERROR: something in calcs_torun which I do not '
                            'recognize: {}'.format(node))
                #self.report('{}{}'.format(type(node[0], node[1])))
                res = None
                continue
            res_all.append(res)
            #print res
            #calc_node = res['output_scf_wc_para'].get_inputs()[0]
            # if run is used, otherwise use labels
            #self.ctx.calcs_res.append(calc_node)
            #self.ctx.calcs_torun.remove(node)
            #print res
        self.ctx.calcs_torun = []
        #return ToContext(last_calc=res)

        '''
        inputs = get_inputs_fleur(code, remote, fleurin, options)
        future = submit(FleurProcess, **inputs)
        self.ctx.loop_count = self.ctx.loop_count + 1
        print 'run FLEUR number: {}'.format(self.ctx.loop_count)
        self.ctx.calcs.append(future)
        '''
        return ToContext(calcs_res=res) #calcs.append(future


    def relaxation_needed(self):
        """
        If the structures should be relaxed, check if their Forces are below a certain
        threshold, otherwise throw them in the relaxation wf.
        """
        self.report('INFO: In relaxation initial_state_CLS workflow (so far nothing to do)')
        if self.ctx.relax:
            # TODO check all forces of calculations
            forces_fine = True
            if forces_fine:
                return True
            else:
                return False
        else:
            return False


    def relax(self):
        """
        Do structural relaxation for certain structures.
        """
        self.report('INFO: In relax initial_state_CLS workflow (so far nothing to do)')
        self.ctx.dos_to_calc = []
        for calc in self.ctx.dos_to_calc:
            pass
            # TODO run relax workflow


    def find_parameters(self):
        """
        If the same parameters shall be used in the calculations you have to
        find some that match. For low error on CLS. therefore use the ones enforced
        or extract from the previous Fleur calculation.
        """
        #self.ctx.ref[elem] = ref_el
        #self.ctx.ref_calcs_torun.append(ref_el)

        # for entry in ref[elem] find parameter node
        for elm, struc in self.ctx.ref.iteritems():
            #print(elm, struc)
            #self.ctx.ref_calcs_torun.append(ref_el)
            pass
            # if parameter node given, extract from there,
            #parameter_dict
            # else
            #extract parameter out of previous calculation
            #parameter_dict = fleurinp.extract_para(element)
            # BE CAREFUL WITH LOs! soc and co


    def run_scfs_ref(self):
        """
        Run SCF-cycles for ref structures, calculations given in certain workflow arrays.
        parameter nodes should be given
        """
        self.report('INFO: In run_scfs_ref initial_state_CLS workflow')

        para = self.ctx.scf_para
        if para == 'default':
            wf_parameter = {}
        else:
            wf_parameter = para
        wf_parameter['serial'] = self.ctx.serial
        # TODO maybe use less resources, or default of one machine
        #wf_parameter['options'] = self.ctx.options 
        wf_parameters = ParameterData(dict=wf_parameter)
        options = ParameterData(dict=self.ctx.options)

        res_all = []
        calcs = {}
        # now in parallel
        #print self.ctx.ref_calcs_torun
        i = 0
        #print(self.ctx.ref_calcs_torun)
        for i, node in enumerate(self.ctx.ref_calcs_torun):
            scf_label = 'cls|scf_wc on ref {}'.format(self.ctx.elements[i])
            scf_description = ('cls|scf of the reference structure of element {}'
                               ''.format(self.ctx.elements[i]))
            #print node
            if isinstance(node, StructureData):
                res = self.submit(fleur_scf_wc, wf_parameters=wf_parameters, structure=node,
                             inpgen = self.inputs.inpgen, fleur=self.inputs.fleur, options=options,
                             label=scf_label, description=scf_description)#
            #elif isinstance(node, FleurinpData):
            #    res = submit(fleur_scf_wc, wf_parameters=wf_parameters, structure=node,
            #                inpgen = self.inputs.inpgen, fleur=self.inputs.fleur)#
            elif isinstance(node, list):#(StructureData, ParameterData)):
                res = self.submit(fleur_scf_wc, wf_parameters=wf_parameters,
                             structure=node[0], calc_parameters=node[1], options=options,
                             inpgen = self.inputs.inpgen, fleur=self.inputs.fleur,
                             label=scf_label, description=scf_description)#
            else:
                self.report('WARNING: something in calcs_torun which I do not reconise: {}'.format(node))
                continue
            label = str('calc_ref{}'.format(i))
            #print(label)
            #calc_node = res['output_scf_wc_para'].get_inputs()[0] # if run is used, otherwise use labels
            self.ctx.ref_labels.append(label)
            calcs[label] = res
            res_all.append(res)
            #print res
            self.ctx.ref_calcs_res.append(res)
            #self.ctx.calcs_torun.remove(node)
            #print res
        self.ctx.ref_calcs_torun = []
        return ToContext(**calcs)

    def handle_scf_failure(self):
        """
        In here we handle all failures from the scf workchain
        """
        '''
        try:
            calculation = self.ctx.calculation
        except Exception:
            self.abort_nowait('the first iteration finished without returning a PwCalculation')
            return

        expected_states = [calc_states.FINISHED, calc_states.FAILED, calc_states.SUBMISSIONFAILED]

        # Done: successful convergence of last calculation
        if calculation.has_finished_ok():
            self.report('converged successfully after {} iterations'.format(self.ctx.iteration))
            self.ctx.restart_calc = calculation
            self.ctx.is_finished = True

        # Abort: exceeded maximum number of retries
        elif self.ctx.iteration >= self.ctx.max_iterations:
            self.report('reached the max number of iterations {}'.format(self.ctx.max_iterations))
            self.abort_nowait('last ran PwCalculation<{}>'.format(calculation.pk))

        # Abort: unexpected state of last calculation
        elif calculation.get_state() not in expected_states:
            self.abort_nowait('unexpected state ({}) of PwCalculation<{}>'.format(
                calculation.get_state(), calculation.pk))

        # Retry: submission failed, try to restart or abort
        elif calculation.get_state() in [calc_states.SUBMISSIONFAILED]:
            self._handle_submission_failure(calculation)

        # Retry: calculation failed, try to salvage or abort
        elif calculation.get_state() in [calc_states.FAILED]:
            self._handle_calculation_failure(calculation)

        # Retry: try to convergence restarting from this calculation
        else:
            self.report('calculation did not converge after {} iterations, restarting'.format(self.ctx.iteration))
            self.ctx.restart_calc = calculation

        '''

        failure = False
        if failure:
            self.ctx.successful = False

        return


    def collect_results(self):
        """
        Collect results from certain calculation, check if everything is fine,
        calculate the wanted quantities. currently all energies are in hartree
        (as provided by Fleur)
        """
        from aiida_fleur.tools.common_fleur_wf_util import determine_formation_energy

        message=('INFO: Collecting results of initial_state_CLS workflow')
        self.report(message)
        # TODO be very careful with core config?
        #from pprint import pprint

        #self.ctx.ref_cl_energies
        all_CLS = {}
        # get results from calc
        calcs = self.ctx.calcs_res
        ref_calcs = self.ctx.ref_calcs_res_given # []
        #print (self.ctx.ref_labels)
        for label in self.ctx.ref_labels:
            calc = self.ctx[label]
            ref_calcs.append(calc)

        # extract_results need the scf workchain calculation node
        total_en, fermi_energies, bandgaps, atomtypes, all_corelevel, log_ref = extract_results([calcs])
        ref_total_en, ref_fermi_energies, ref_bandgaps, ref_atomtypes, ref_all_corelevel, log = extract_results(ref_calcs)

        if log_ref or log:
            self.ctx.successful = False
            self.ctx.warnings.append(log_ref)
            self.ctx.warnings.append(log)
        #print(all_corelevel)
        #print(ref_all_corelevel)

        ref_cl_energies = {}
        cl_energies = {}

        #first substract efermi from corelevel of reference structures
        # TODO check if both values, corelevel and efermi are in eV
        for compound, atomtypes_list in ref_atomtypes.iteritems():
            # atomtype_list contains a list of dicts of all atomtypes from compound x
            # get corelevels of compound x
            cls_all_atomtyps = ref_all_corelevel[compound]
            for i, atomtype in enumerate(atomtypes_list):
                #atomtype a dict which contains one atomtype
                elm = atomtype.get('element', None)
                cls_atomtype = cls_all_atomtyps[i][0]
                ref_cl_energies[elm] = []
                ref_cls = []
                for corelevel in cls_atomtype['corestates']:
                    ref_cls.append(corelevel['energy'] - ref_fermi_energies[compound])
                ref_cl_energies[elm].append(ref_cls)

        #print('ref_cl energies')
        #print(ref_cl_energies)
        #print(all_corelevel)

        #now substract efermi from corelevel of compound structure
        #and calculate core level shifts
        for compound, cls_atomtypes_list in all_corelevel.iteritems():
            #init, otherwise other types will override
            for i, atomtype in enumerate(atomtypes[compound]):
                elm = atomtype.get('element', None)
                cl_energies[elm] = []
                all_CLS[elm] = []

            #now fill
            for i, atomtype in enumerate(atomtypes[compound]):
                elm = atomtype.get('element', None)
                #print elm
                cls_atomtype = cls_atomtypes_list[i]
                corelevels = []
                for corelevel in cls_atomtype[0]['corestates']:
                    correct_cl = corelevel['energy']-fermi_energies[compound]
                    corelevels.append(correct_cl)
                cl_energies[elm].append(corelevels)

                #now calculate CLS
                ref = ref_cl_energies.get(elm,[0])[-1]# We just use one (last) atomtype
                #of elemental reference (in general might be more complex,
                #since certain elemental unit cells could have several atom types (graphene))
                corelevel_shifts = []
                #TODO shall we store just one core-level shift per atomtype?
                for i, corelevel in enumerate(cl_energies[elm][-1]):
                    corelevel_shifts.append(corelevel - float(ref[i]))
                all_CLS[elm].append(corelevel_shifts)

        # calculate formation energy
        #determine_formation_energy({'BeW' : 2, 'Be2W' : 2.5}, {'Be' : 1, 'W' : 1})
        # to normalize ref,
        # from.split(012345678910)
        # devide total energy by number of atoms

        # Formation energy calculation is ony possible if all elementals of the structure
        # have been calculated.
        # convert total_en dict to list, why?
        total_en_list =[]
        for key, val in total_en.iteritems():
            total_en_list.append([key, val])
        if self.ctx.calculate_formation_energy:
            # the reference total energy is for the whole structure with several atoms,
            # we need it per atom
            ref_total_en_norm = {}
            for key, val in ref_total_en.iteritems():
                elm_dict = get_natoms_element(key)
                ref_total_en_norm[elm_dict.keys()[0]] = 1.0* val/elm_dict.values()[0]
            #print ref_total_en_norm
            #print total_en
            
            formation_energy, form_dict = determine_formation_energy(total_en, ref_total_en_norm)
        else:
            formation_energy = [[]]

        # TODO make simpler format of atomtypes for node
        # TODO write corelevel explanation/coresetup in a format like 4f7/2
        #TODO ? also get total energies?
        return cl_energies, all_CLS, ref_cl_energies, fermi_energies, bandgaps, ref_fermi_energies, ref_bandgaps, atomtypes, ref_atomtypes, formation_energy[0], total_en, ref_total_en

    def return_results(self):
        """
        return the results of the calculations
        """
        # TODO more output, info here

        #print corelevel shifts were calculated bla bla
        cl, cls, ref_cl, efermi, gap, ref_efermi, ref_gap, at, at_ref, formE, tE, tE_ref = self.collect_results()

        if self.ctx.errors:
            self.ctx.warnings.append(self.ctx.errors)

        outputnode_dict = {}

        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['workflow_version'] = self._workflowversion
        outputnode_dict['warnings'] = self.ctx.warnings
        outputnode_dict['successful'] = self.ctx.successful
        outputnode_dict['material'] = efermi.keys()[0]
        outputnode_dict['corelevel_energies'] = cl
        outputnode_dict['corelevel_energies_units'] = 'htr'#'eV'
        outputnode_dict['reference_corelevel_energies'] = ref_cl
        outputnode_dict['reference_corelevel_energies_units'] = 'htr'#'eV'
        outputnode_dict['reference_fermi_energy'] = ref_efermi.values()
        outputnode_dict['reference_fermi_energy_des'] = ref_efermi.keys()
        outputnode_dict['fermi_energy'] = efermi.values()[0]
        outputnode_dict['fermi_energy_units'] = 'htr'
        outputnode_dict['corelevelshifts'] = cls
        outputnode_dict['corelevelshifts_units'] = 'htr'
        outputnode_dict['binding_energy_convention'] = 'negativ'
        #outputnode_dict['coresetup'] = []#cls
        #outputnode_dict['reference_coresetup'] = []#cls
        outputnode_dict['bandgap'] = gap.values()[0]
        outputnode_dict['bandgap_units'] = 'htr'
        outputnode_dict['reference_bandgaps'] = ref_gap.values()
        outputnode_dict['reference_bandgaps_des'] = ref_gap.keys()
        outputnode_dict['atomtypes'] = at
        outputnode_dict['formation_energy'] = formE
        outputnode_dict['formation_energy_units'] = 'eV/atom'
        outputnode_dict['total_energy'] = tE.values()[0]
        outputnode_dict['total_energy_units'] = 'eV'
        outputnode_dict['total_energy_ref'] = tE_ref.values()
        outputnode_dict['total_energy_ref_des'] = tE_ref.keys()
        #outputnode = ParameterData(dict=outputnode_dict)

        # To have to ouput node linked to the calculation output nodes
        outnodedict = {}
        outnode = ParameterData(dict=outputnode_dict)
        outnodedict['results_node'] = outnode

        # TODO: bad design, put in workfunction and make bullet proof.
        calc = self.ctx.calcs_res
        calc_dict = calc.get_outputs_dict()['output_scf_wc_para']
        outnodedict['input_structure']  = calc_dict

        for label in self.ctx.ref_labels:
            calc = self.ctx[label]
            calc_dict = calc.get_outputs_dict()['output_scf_wc_para']
            outnodedict[label]  = calc_dict

        outdict = create_initcls_result_node(**outnodedict)

        #outdict = {}
        #outdict['output_initial_cls_wc_para'] = outputnode
        #print outdict
        for k, v in outdict.iteritems():
            self.out(k, v)
        msg = ('INFO: Initial_state_CLS workflow Done')
        self.report(msg)

    def control_end_wc(self, errormsg):
        """
        Controled way to shutdown the workchain. will initalize the output nodes
        The shutdown of the workchain will has to be done afterwards
        """
        self.ctx.successful = False
        self.ctx.abort = True
        self.report(errormsg) # because return_results still fails somewhen
        self.ctx.errors.append(errormsg)
        self.return_results()
        
        return        
       


@wf
def create_initcls_result_node(**kwargs):
    """
    This is a pseudo wf, to create the rigth graph structure of AiiDA.
    This wokfunction will create the output node in the database.
    It also connects the output_node to all nodes the information commes from.
    So far it is just also parsed in as argument, because so far we are to lazy
    to put most of the code overworked from return_results in here.
    """

    outdict = {}
    outpara = kwargs.get('results_node', {})
    outdict['output_initial_cls_wc_para'] = outpara.clone()
    # copy, because we rather produce the same node twice then have a circle in the database for now...
    #output_para = args[0]
    #return {'output_eos_wc_para'}
    return outdict



def querry_for_ref_structure(element_string):
    """
    This methods finds StructureData nodes with the following extras:
    extra.type = 'bulk', # Should be done by looking at pbc, but I could not
    get querry to work.
    extra.specific = 'reference',
    'extra.elemental' = True,
    extra.structure = element_string

    param: element_string: string of an element
    return: the latest StructureData node that was found
    """

    #query db
    q = QueryBuilder()
    q.append(StructureData,
             filters = {'extras.type' : {'==' : 'bulk'},
                        'extras.specification' : {'==' : 'reference'},
                        'extras.elemental' : {'==' : True},
                        'extras.element' : {'==' : element_string}
                        })
    q.order_by({StructureData : 'ctime'})#always use the most recent
    structures = q.all()

    if structures:
        return structures[-1][0]
    else:
        return None


def fleur_calc_get_structure(calc_node):
    """
    Get the AiiDA data structure from a fleur calculations
    """
    #get fleurinp
    fleurinp = calc_node.inp.fleurinpdata
    structure = fleurinp.get_structuredata(fleurinp)
    return structure

def extract_results(calcs):
    """
    Collect results from certain calculation, check if everything is fine,
    calculate the wanted quantities.

    params: calcs : list of scf workchains nodes
    """

    from aiida_fleur.tools.extract_corelevels import extract_corelevels
    log = []
    calc_uuids = []
    for calc in calcs:
        #print(calc)
        try:
            calc_uuids.append(calc.get_outputs_dict()['output_scf_wc_para'].get_dict()['last_calc_uuid'])
        except: #TODO what error
            logmsg = ('ERROR: No output_scf_wc_para node found for "last_calc_uuid"'
                      'key in it for calculation: {}'.format(calc))
            log.append(logmsg)
            continue
        #calc_uuids.append(calc['output_scf_wc_para'].get_dict()['last_calc_uuid'])
    #print(calc_uuids)

    all_corelevels = {}
    fermi_energies = {}
    bandgaps = {}
    all_atomtypes = {}
    total_energy = {}

    # more structures way: divide into this calc and reference calcs.
    # currently the order in calcs is given, but this might change if you submit
    # check if calculation pks belong to successful fleur calculations
    for uuid in calc_uuids:
        calc = load_node(uuid)
        if not isinstance(calc, FleurCalc):
            #raise ValueError("Calculation with pk {} must be a FleurCalculation".format(pk))
            # log and continue
            continue
        if calc.get_state() == 'FINISHED':
            # get out.xml file of calculation
            outxml = calc.out.retrieved.folder.get_abs_path('path/out.xml')
            #print outxml
            corelevels, atomtypes = extract_corelevels(outxml)
            #all_corelevels.append(core)
            #print('corelevels: {}'.format(corelevels))
            #print('atomtypes: {}'.format(atomtypes))
            #for i in range(0,len(corelevels[0][0]['corestates'])):
            #    print corelevels[0][0]['corestates'][i]['energy']

            #TODO how to store?
            efermi = calc.res.fermi_energy
            #print efermi
            bandgap = calc.res.bandgap
            te = calc.res.energy
            #total_energy = calc.res.total_energy
            #total_energy_units = calc.res.total_energy_units

            # TODO: maybe different, because it is prob know from before
            #fleurinp = calc.inp.fleurinpdata
            #structure = fleurinp.get_structuredata(fleurinp)
            #compound = structure.get_formula()
            #print compound
            #number = '{}'.format(i)
            #fermi_energies[number] = efermi
            #bandgaps[number] = bandgap
            #all_atomtypes[number] = atomtypes
            #all_corelevels[number] = corelevels
            #all_total_energies[number] = total_energy
        else:
            # log and continue
            te = float('nan')
            bandgap = float('nan')
            efermi = float('nan')
            corelevels = [float('nan')]
            atomtypes = [float('nan')]
            logmsg = "ERROR: Fleur Calculation with uuid {} was not in in state FINISHED".format(uuid)
            log.append(logmsg)
            #continue
            #raise ValueError("Calculation with pk {} must be in state FINISHED".format(pk))

        # TODO: maybe different, because it is prob know from before
        fleurinp = calc.inp.fleurinpdata
        structure = fleurinp.get_structuredata_nwf()
        compound = structure.get_formula()
        #print compound
        fermi_energies[compound] = efermi
        bandgaps[compound] = bandgap
        all_atomtypes[compound] = atomtypes
        all_corelevels[compound] = corelevels
        total_energy[compound] = te
        #fermi_energies = efermi
        #bandgaps = bandgap
        #all_atomtypes = atomtypes
        #all_corelevels = corelevels

    return total_energy, fermi_energies, bandgaps, all_atomtypes, all_corelevels, log
    #TODO validate results and give some warnings

    # check bandgaps, if not all metals, throw warnings:
    # bandgap and efermi prob wrong, which makes some results meaningless

    # check fermi energy differences, correct results for fermi energy diff
    # ggf TODO make a raw core-level and core-level to fermi energy variable
    #TODO to what reference energy? or better not to fermi, but first unocc? (add bandgap)

    #Style: {atomtype : listof all corelevel, atomtype_coresetup... }
    #ie: { 'W-1' : [shift_1s, ... shift 7/2 4f],
    #      'W-1_coreconfig' : ['1s','2s',...],
    #      'W-2' : [...], 'Be-1': [], ...} #all in eV!
    #self.ctx.CLS = {}
    #self.ctx.cl_energies = {}# same style as CLS only energy <-> shift

    #Style: {'Compound' : energy, 'ref_x' : energy , ...}
    #i.e {'Be12W' : 0.0, 'Be' : 0.104*htr_eV , 'W' : 0.12*htr_eV} # all in eV!
    #self.ctx.fermi_energies = {}


def get_ref_from_group(element, group):
    """
    Return a structure data node from a given group for a given element.
    (quit creedy, done straighforward)

    params: group: group name or pk
    params: element: string with the element i.e 'Si'

    returns: AiiDA StructureData node
    """

    report = []

    try:
        group_pk = int(group)
    except ValueError:
        group_pk = None
        group_name = group

    if group_pk is not None:
        try:
            str_group = Group(dbgroup=group_pk)
        except NotExistent:
            str_group = None
            message = ('You have to provide a valid pk for a Group of'
                       'structures or a Group name. Reference key: "group".'
                       'given pk= {} is not a valid group'
                       '(or is your group name integer?)'.format(group_pk))
            report.append(message)
    else:
        try:
            str_group = Group.get_from_string(group_name)
        except NotExistent:
            str_group = None
            message = ('You have to provide a valid pk for a Group of'
                       'structures or a Group name. Wf_para key: "struc_group".'
                       'given group name= {} is not a valid group'
                       '(or is your group name integer?)'.format(group_name))
            report.append(message)
            #abort_nowait('I abort, because I have no structures to calculate ...')

    stru_nodes = str_group.nodes
    #n_stru = len(stru_nodes)

    structure = None

    for struc in stru_nodes:
        formula = struc.get_formula()
        eformula = formula.translate(None, digits) # remove digits, !python3 differs
        if eformula == element:
            return struc, report

    report.append('Structure node for element {} not found in group {}'
                  ''.format(element, group))

    return structure, report


def get_para_from_group(element, group):
    """
    get structure node for a given element from a given group of structures
    (quit creedy, done straighforward)
    """

    report = []

    try:
        group_pk = int(group)
    except ValueError:
        group_pk = None
        group_name = group

    if group_pk is not None:
        try:
            para_group = Group(dbgroup=group_pk)
        except NotExistent:
            para_group = None
            message = ('You have to provide a valid pk for a Group of '
                       'parameters or a Group name. Reference key: "group".'
                       'given pk= {} is not a valid group'
                       '(or is your group name integer?)'.format(group_pk))
            report.append(message)
    else:
        try:
            para_group = Group.get_from_string(group_name)
        except NotExistent:
            para_group = None
            message = ('You have to provide a valid pk for a Group of '
                       'parameters or a Group name. Wf_para key: "para_group".'
                       'given group name= {} is not a valid group'
                       '(or is your group name integer?)'.format(group_name))
            report.append(message)
            #abort_nowait('I abort, because I have no structures to calculate ...')

    para_nodes = para_group.nodes
    #n_stru = len(para_nodes)

    parameter = None

    for para in para_nodes:
        formula = para.get_extras().get('element', None)
        #eformula = formula.translate(None, digits) # remove digits, !python3 differs
        if formula == element:
            return para, report

    report.append('Parameter node for element {} not found in group {}'
                  ''.format(element, group))

    return parameter, report



def clshifts_to_be(coreleveldict, reference_dict):
    """
    This methods converts corelevel shifts to binding energies,
    if a reference is given.
    These cann than be used for plotting.

    i.e:

    reference = {'W' : {'4f7/2' : [124],
                     '4f5/2' : [102]},
              'Be' : {'1s': [117]}}
    corelevels = {'W' : {'4f7/2' : [0.4, 0.3, 0.4 ,0.1],
                     '4f5/2' : [0, 0.3, 0.4, 0.1]},
              'Be' : {'1s': [0, 0.2, 0.4, 0.1, 0.3]}}
    """

    return_corelevel_dict = {}
    for elem, corelevel_dict in coreleveldict.iteritems():
        ref_el = reference_dict.get(elem, {})
        return_corelevel_dict[elem] = {}
        for corelevel_name, corelevel_list in corelevel_dict.iteritems():
            ref_cl = ref_el.get(corelevel_name, [])
            be_all = []
            nref = len(ref_cl)
            ncl = len(corelevel_list)
            if nref == ncl:
                for i, corelevel in enumerate(corelevel_list):
                    be = corelevel + ref_cl[i]
                    be_all.append(be)
            else:
                for corelevel in corelevel_list:
                    be = corelevel + ref_cl[0]
                    be_all.append(be)
            return_corelevel_dict[elem][corelevel_name] = be_all

    return return_corelevel_dict

