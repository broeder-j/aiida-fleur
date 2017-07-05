#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This is the worklfow 'corelevel' using the Fleur code, which calculates Binding
energies and corelevel shifts with different methods.
'divide and conquer'
"""
#TODO parsing of eigenvalues of LOS!
#TODO error handling of scf
#TODO USE SAME PARAMETERS! (maybe extract method for fleurinp needed)

from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()
    
#import os.path
from aiida.orm import Code, DataFactory, CalculationFactory, load_node
from aiida.work.workchain import WorkChain
from aiida.work.run import submit
from aiida.work.run import async as asy
from aiida.work.workchain import ToContext
from aiida.work.process_registry import ProcessRegistry
from aiida_fleur.workflows.scf import fleur_scf_wc

from aiida_fleur.calculation.fleur import FleurCalculation
from aiida_fleur.data.fleurinpmodifier import FleurinpModifier
from aiida.work.workchain import  if_ #while_,
from aiida_fleur_ad.util.extract_corelevels import extract_corelevels
from aiida.common.exceptions import NotExistent

StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
RemoteData = DataFactory('remote')
FleurinpData = DataFactory('fleur.fleurinp')
FleurProcess = FleurCalculation.process()
FleurCalc = CalculationFactory('fleur.fleur')

htr_to_eV = 1

class fleur_initial_cls_wc(WorkChain):
    '''
    Turn key solution for the calculation of core level shift and Binding energies
    
    '''
    # wf_Parameters: ParameterData, 
    '''
    'method' : ['initial', 'full_valence ch', 'half_valence_ch', 'ch', ...]
    'Bes' : [W4f, Be1s]
    'CLS' : [W4f, Be1s]
 toms' : ['all', 'postions' : []]
    #'references' : ['calculate', and use # calculate : 'all' , or 'calculate' : ['W', 'Be']
    'references' : { 'W': [calc/ouputnode or  fleurinp, or structure data or structure data + Parameter  ], 'Be' : }
    'scf_para' : {...}, 'default' 
    'relax' : True
    'relax_mode': ['Fleur', 'QE Fleur', 'QE']
    'relax_para' : {...}, 'default' 
    'calculate_doses' : False
    'dos_para' : {...}, 'default' 
    '''
    '''
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
    '''
    
    _workflowversion = "0.0.1"
    _default_wf_para = {#'references' : {'calculate' : 'all'},
                        'structure_ref' : {},
                        'relax' : True,
                        'relax_mode': 'QE Fleur',
                        'relax_para' : 'default',
                        'scf_para' : 'default',
                        'same_para' : True,
                        'resources' : {"num_machines": 1},
                        'walltime_sec' : 10*30,
                        'queue_name' : None,
                        'serial' : True}    
    def __init__(self, *args, **kwargs):
        super(fleur_initial_cls_wc, self).__init__(*args, **kwargs)    
    '''
    def get_defaut_wf_para(self):
        return self._default_wf_para
     '''     
    @classmethod
    def define(cls, spec):
        super(fleur_initial_cls_wc, cls).define(spec)
        spec.input("wf_parameters", valid_type=ParameterData, required=False,
                   default=ParameterData(dict={ 
                        #'references' : {'calculate' : 'all'}, 
                        'references' : {},
                        'relax' : True,
                        'relax_mode': 'QE Fleur',
                        'relax_para' : 'default',
                        'scf_para' : 'default',
                        'same_para' : True,
                        'resources' : {"num_machines": 1},
                        'walltime_sec' : 10*60,
                        'queue_name' : None,
                        'serial' : True}))#TODO_default_wf_para out of here#
        spec.input("fleurinp", valid_type=FleurinpData, required=False)
        spec.input("fleur", valid_type=Code, required=True)
        spec.input("inpgen", valid_type=Code, required=False)        
        spec.input("structure", valid_type=StructureData, required=False)
        spec.input("calc_parameters", valid_type=ParameterData, required=False)
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
        spec.dynamic_output()
        #spec.dynamic_input()

    def check_input(self):
        """
        Init same context and check what input is given if it makes sence
        """
        ### input check ### ? or done automaticly, how optional?

        msg=("INFO: Started inital_state_CLS workflow version {} "
             "Workchain node identifiers: {}"
              "".format(self._workflowversion, ProcessRegistry().current_calc_node))
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
        self.ctx.successful = False
        self.ctx.warnings = []
        self.ctx.errors = []
        self.ctx.ref = {}
        
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
        self.ctx.resources = wf_dict.get('resources', default.get('resources'))
        self.ctx.walltime_sec = wf_dict.get('walltime_sec', default.get('walltime_sec'))
        self.ctx.queue = wf_dict.get('queue_name', default.get('queue_name'))
        # check if inputs given make sense
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
            if not 'inpgen' in inputs:
                error = 'ERROR: StructureData was provided, but no inpgen code was provided'
                self.ctx.errors.append(error)
                self.abort_nowait(error)
            if 'calc_parameters' in inputs:
                self.ctx.calcs_torun.append((inputs.get('calc_parameters'), inputs.get('structure')))
                #print('here2')
            else:
                self.ctx.calcs_torun.append(inputs.get('structure'))
                #print('here3')
        else:
            error = 'ERROR: No StructureData nor FleurinpData was provided'
            #print(error)
            self.ctx.errors.append(error)
            self.abort_nowait(error)
        self.report('INFO: elements in structure: {}'.format(self.ctx.elements))
        
        
    def get_references(self):
        """
        To calculate a CLS in inital state approx, we need reference calculations
        to the Elemental crystals. First it is checked if the user has provided them
        Second the database is checked, if there are structures with certain extras.
        Third the COD database is searched for the elemental Cystal structures.
        If some referneces are not found stop here.
        Are there already calculation of these 'references', ggf use them.
        We do not put these calculation in the calculation queue yet because we
        need specific parameters for them
        """

        self.report('INFO: In Get_references inital_state_CLS workflow')   

        #references = self.inputs.wf_parameters.get_dict().get('references', {'calculate' : 'all'})
        references = self.inputs.wf_parameters.get_dict().get('references', {})
        # should be of the form of
        #'references' : { 'W': calc, outputnode of workflow or fleurinp, 
                         #or structure data or (structure data + Parameter), 
        #                 'Be' : ...}
        
        self.ctx.ref = {}
        self.ctx.abort = False

        struc_group = references.get('group', None)
        para_group = references.get('para_group', None)

        #TODO better checks if ref makes sense?
        
        # get specific element reference if given override
        for elem in self.ctx.elements:
            #to_calc[elem] = 'find' 
            ref_el = references.get(elem, None)
            if ref_el:
                try:
                    ref_el_node = load_node(ref_el)
                except:
                    ref_el_node = None
                    self.report('ERROR: The reference node (id or uuid) provided: {} for '
                                'element: {} could not be loaded with load_node'.format(ref_el, elem))
                    self.ctx.abort = True
                if isinstance(ref_el_node, (StructureData, ParameterData)):
                    #self.ctx.ref[elem] = ref_el
                    #enforced parameters, add directly to run queue
                    self.ctx.ref_calcs_torun.append(ref_el_node)
                elif isinstance(ref_el_node, FleurCalc):
                    #extract from fleur calc TODO
                    self.ctx.ref_cl_energies[elem] = {}
                elif isinstance(ref_el_node, ParameterData):
                    #extract from workflow output TODO
                    self.ctx.ref_cl_energies[elem] = {}             
                elif isinstance(ref_el_node, FleurinpData):
                    # add to calculations
                    #enforced parameters, add directly to run queue
                    self.ctx.ref_calcs_torun.append(ref_el_node)
                    #self.ctx.ref[elem] = ref_el
                elif isinstance(ref_el_node, StructureData):
                    self.ctx.ref[elem] = ref_el_node
                #elif isinstance(ref_el, initial_state_CLS):
                #    extract TODO
                else:
                    error = ("ERROR: I do not know what to do with this given reference"
                             "{} for element {}".format(ref_el, elem))
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
                    self.ctx.ref_calcs_torun.append((structure, parameter))
                elif structure:
                    self.ctx.ref[elem] = structure
                    self.ctx.ref_calcs_torun.append(structure)
                else:
                    pass # report not found?
            else: # no ref given, we have to look for it.
                structure = querry_for_ref_structure(elem)
                if structure:
                    self.ctx.ref[elem] = structure
                    self.ctx.ref_calcs_torun.append(structure)# tempoary later check parameters
                else: #not found
                    error = ("ERROR: Reference structure for element: {} not found."
                             "checkout the 'querry_for_ref_structure' method."
                             "to see what extras are querried for.".format(elem))
                    #print(error)
                    self.ctx.errors.append(error)
                    self.ctx.abort = True
                    self.report(error)
        if self.ctx.abort:
            error = ('ERROR: Something was wrong with the reference input provided. '
                    'I cannot calculate from the input, or what I have found '
                    'what you want me to do. Please check the workchain report'
                    'for details.')
            self.abort_nowait(error)

        #print('self.ctx.ref: {} '.format(self.ctx.ref))
        #StructureData 
        #ParameterData
        #FleurinpData
        #FleurCalc
        
        # check if a structureData for these elements was given
        #if yes add to ref_calc to run
        #was also a prameter node given for the element?
        #yes run with these
        #no was on given for the host structure, extract element parameternode
        
        #else use parameters extracted from host calculation # TODO
        
        #check if there is a structure from this element in the database with extras: 
        # with extra.type = 'bulk', extra.specific = 'reference', 'extra.elemental' = True, extra.structure = 'W'
        # check if input parameter node values for the calculation are the same.
        
        #if yes, if a calculation exists use that result
        #else do a calculation on that structure as above
 

    
    def run_fleur_scfs(self):
        """
        Run SCF-cycles for all structures, calculations given in certain workflow arrays.
        """
        self.report('INFO: In run_fleur_scfs inital_state_CLS workflow')        
        #from aiida.work import run, async, 
        #TODO if submiting of workdlows work, use that. 
        #or run them with async (if youy know how to extract results) 
        
        para = self.ctx.scf_para
        if para == 'default': 
            wf_parameter = {}
        else:
            wf_parameter = para
        
        wf_parameter['queue_name'] = self.ctx.queue
        wf_parameter['serial'] = self.ctx.serial
        wf_parameters =  ParameterData(dict=wf_parameter)
        res_all = []
        # for each calulation in self.ctx.calcs_torun #TODO what about wf params?
        #print self.ctx.calcs_torun
        for node in self.ctx.calcs_torun:
            #print node
            if isinstance(node, StructureData):
                res = fleur_scf_wc.run(wf_parameters=wf_parameters, structure=node, 
                            inpgen = self.inputs.inpgen, fleur=self.inputs.fleur)#
            elif isinstance(node, FleurinpData):
                res = fleur_scf_wc.run(wf_parameters=wf_parameters, structure=node, 
                            inpgen = self.inputs.inpgen, fleur=self.inputs.fleur)#
            elif isinstance(node, (StructureData, ParameterData)):
                res = fleur_scf_wc.run(wf_parameters=wf_parameters, calc_parameters=node(1), structure=node(0), 
                            inpgen = self.inputs.inpgen, fleur=self.inputs.fleur)#
            else:
                self.report('ERROR: something in calcs_torun which I do not recognize: {}'.format(node))
                self.report('{}{}'.format(type(node[0], node[1])))
                continue
            res_all.append(res)
            #print res
            calc_node = res['output_scf_wc_para'].get_inputs()[0] # if run is used, otherwise use labels
            self.ctx.calcs_res.append(calc_node)
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
        #return ToContext(last_calc=res) #calcs.append(future

        
    def relaxation_needed(self):
        """
        If the structures should be relaxed, check if their Forces are below a certain 
        threshold, otherwise throw them in the relaxation wf.
        """
        self.report('INFO: In relaxation inital_state_CLS workflow (so far nothing to do)')
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
        self.report('INFO: In relax inital_state_CLS workflow (so far nothing to do)')  
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
        self.report('INFO: In run_scfs_ref inital_state_CLS workflow')        
        #from aiida.work import run, async, 
        #TODO if submiting of workdlows work, use that. 
        #or run them with async (if youy know how to extract results) 
        
        para = self.ctx.scf_para
        if para == 'default': 
            wf_parameter = {}
        else:
            wf_parameter = para
        wf_parameter['serial'] = self.ctx.serial
        wf_parameter['queue_name'] = self.ctx.queue
        wf_parameters =  ParameterData(dict=wf_parameter)
        res_all = []
        calcs = {}
        # now in parallel
        #print self.ctx.ref_calcs_torun
        i = 0
        for node in self.ctx.ref_calcs_torun:
            #print node
            i = i+1
            if isinstance(node, StructureData):
                res = asy(fleur_scf_wc, wf_parameters=wf_parameters, structure=node, 
                            inpgen = self.inputs.inpgen, fleur=self.inputs.fleur)#
            elif isinstance(node, FleurinpData):
                res = asy(fleur_scf_wc, wf_parameters=wf_parameters, structure=node, 
                            inpgen = self.inputs.inpgen, fleur=self.inputs.fleur)#
            elif isinstance(node, (StructureData, ParameterData)):
                res = asy(fleur_scf_wc, wf_parameters=wf_parameters, calc_parameters=node(1), structure=node(0), 
                            inpgen = self.inputs.inpgen, fleur=self.inputs.fleur)#
            else:
                print('something in calcs_torun which I do not reconise: {}'.format(node))
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
      
        '''
        # for each calulation in self.ctx.calcs_torun #TODO what about wf params?
        print self.ctx.ref_calcs_torun
        for node in self.ctx.ref_calcs_torun:
            print node
            if isinstance(node, StructureData):
                res = fleur_scf_wc.run(wf_parameters=wf_parameters, structure=node, 
                            inpgen = self.inputs.inpgen, fleur=self.inputs.fleur)#
            elif isinstance(node, FleurinpData):
                res = fleur_scf_wc.run(wf_parameters=wf_parameters, structure=node, 
                            inpgen = self.inputs.inpgen, fleur=self.inputs.fleur)#
            elif isinstance(node, (StructureData, ParameterData)):
                res = fleur_scf_wc.run(wf_parameters=wf_parameters, calc_parameters=node(1), structure=node(0), 
                            inpgen = self.inputs.inpgen, fleur=self.inputs.fleur)#
            else:
                print('something in calcs_torun which I do not reconise: {}'.format(node))
                continue
            res_all.append(res)
            print res  
            self.ctx.ref_calcs_res.append(res)
            #self.ctx.calcs_torun.remove(node)
            #print res    
        self.ctx.ref_calcs_torun = []
        #return ToContext(last_calc=res)
        '''
    def handle_scf_failure(self):
        """
        In here we handle all failures from the scf workchain
        """
        pass
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
        return        

        
    def collect_results(self):
        """
        Collect results from certain calculation, check if everything is fine, 
        calculate the wanted quantities. currently all energies are in hartree (as provided by Fleur)
        """
        message=('INFO: Collecting results of inital_state_CLS workflow')
        self.report(message)        
        # TODO be very careful with core config?
        #from pprint import pprint
        
        #self.ctx.ref_cl_energies
        all_CLS = {}
        # get results from calc
        calcs = self.ctx.calcs_res
        ref_calcs = []#self.ctx.ref_calcs_res
        #print (self.ctx.ref_labels)
        for label in self.ctx.ref_labels:
            calc = self.ctx[label]
            ref_calcs.append(calc)
        
        #print('ref_calcs')
        #print ref_calcs
        #print('calcs')
        #print calcs
        # extract_results need the scf workchain calculation node
        fermi_energies, bandgaps, atomtypes, all_corelevel = extract_results(calcs)
        ref_fermi_energies, ref_bandgaps, ref_atomtypes, ref_all_corelevel = extract_results(ref_calcs)

        #print(all_corelevel)
        #print(ref_all_corelevel)
        
        ref_cl_energies = {}
        cl_energies = {}
        
        #first substract efermi from corelevel of reference structures
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
                    ref_cls.append(corelevel['energy']-ref_fermi_energies[compound])
                ref_cl_energies[elm].append(ref_cls)
        
        #print('ref_cl energies')
        #print(ref_cl_energies)
        #pprint(all_corelevel)
        
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
                ref = ref_cl_energies[elm][-1]# We just use one (last) atomtype
                #of elemental reference (in general might be more complex,
                #since certain elemental unit cells could have several atom types (graphene))
                corelevel_shifts = []
                #TODO shall we store just one core-level shift per atomtype?
                for i, corelevel in enumerate(cl_energies[elm][-1]):
                    corelevel_shifts.append(corelevel - float(ref[i]))
                all_CLS[elm].append(corelevel_shifts)
        
        # TODO make simpler format of atomtypes for node
        # TODO write corelevel explanation/coresetup in a format like 4f7/2 
        #TODO ? also get total energies?
        return cl_energies, all_CLS, ref_cl_energies, fermi_energies, bandgaps, ref_fermi_energies, ref_bandgaps, atomtypes, ref_atomtypes
        
    def return_results(self):
        """
        return the results of the calculations
        """
        # TODO more output, info here
        
        #print corelevel shifts were calculated bla bla
        cl, cls, ref_cl, efermi, gap, ref_efermi, ref_gap, at, at_ref =  self.collect_results()
        
        outputnode_dict ={}
        
        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['warnings'] = self.ctx.warnings               
        outputnode_dict['successful'] = self.ctx.successful
        outputnode_dict['corelevel_energies'] = cl #self.ctx.cl_energies
        outputnode_dict['reference_corelevel_energies'] = ref_cl #self.ctx.cl_energies
        outputnode_dict['fermi_energy'] = efermi #self.ctx.fermi_energies               
        outputnode_dict['corelevelshifts'] = cls #self.ctx.CLS
        outputnode_dict['coresetup'] = []#cls
        outputnode_dict['reference_coresetup'] = []#cls
        outputnode_dict['bandgap'] = gap#self.ctx.bandgaps
        outputnode_dict['reference_bandgaps'] = ref_gap#self.ctx.bandgaps
        outputnode_dict['atomtypes'] = at#self.ctx.atomtypes
        #print outputnode_dict
        outputnode = ParameterData(dict=outputnode_dict)
        outdict = {}
        outdict['output_inital_cls_wc_para'] = outputnode
        #print outdict
        for k, v in outdict.iteritems():
            self.out(k, v)
        msg=('INFO: Inital_state_CLS workflow Done')
        self.report(msg)

    def create_new_fleurinp(self):
        """
        create a new fleurinp from the old with certain parameters
        """
        # TODO allow change of kpoint mesh?, tria?
        wf_dict = self.inputs.wf_parameters.get_dict()
        nkpts = wf_dict.get('nkpts', 500) 
        # how can the user say he want to use the given kpoint mesh, ZZ nkpts : False/0
        sigma = wf_dict.get('sigma', 0.005)
        emin = wf_dict.get('emin', -0.30)
        emax = wf_dict.get('emax', 0.80)
      
        fleurmode = FleurinpModifier(self.inputs.fleurinp)

        #change_dict = {'band': True, 'ndir' : -1, 'minEnergy' : self.inputs.wf_parameters.get_dict().get('minEnergy', -0.30000000), 
        #'maxEnergy' :  self.inputs.wf_parameters.get_dict().get('manEnergy','0.80000000'), 
        #'sigma' :  self.inputs.wf_parameters.get_dict().get('sigma', '0.00500000')}
        change_dict = {'band': True, 'ndir' : 0, 'minEnergy' : emin,
                       'maxEnergy' : emax, 'sigma' : sigma} #'ndir' : 1, 'pot8' : True
        
        fleurmode.set_inpchanges(change_dict)

        if nkpts:
            fleurmode.set_nkpts(count=nkpts)
            #fleurinp_new.replace_tag()
        
        fleurmode.show(validate=True, display=False) # needed?
        fleurinp_new = fleurmode.freeze()
        self.ctx.fleurinp1 = fleurinp_new
        #print(fleurinp_new)
        #print(fleurinp_new.folder.get_subfolder('path').get_abs_path(''))

        
    def run_fleur(self):
        '''
        run a fleur calculation
        '''
        FleurProcess = FleurCalculation.process()
        inputs = {}
        inputs = self.get_inputs_fleur()
        #print inputs
        future = submit(FleurProcess, **inputs)
        #print 'run Fleur in band workflow'

        return ToContext(last_calc=future)

def querry_for_ref_structure(element_string):
    """
    This methods finds StructureData nodes with the following extras:
    extra.type = 'bulk', # Should be done by looking at pbc, but I could not get querry to work.
    extra.specific = 'reference', 
    'extra.elemental' = True, 
    extra.structure = element_string
    
    param: element_string: string of an element
    return: the latest StructureData node that was found
    
    """
    from aiida.orm.querybuilder import QueryBuilder

    #query db
    q = QueryBuilder()
    q.append(StructureData,
        filters = {
            'extras.type' : {'==' : 'bulk'},
            'extras.specification' : {'==' : 'reference'},
            'extras.elemental' : {'==' : True},
            'extras.element' : {'==' : element_string}
            }
        )
    q.order_by({StructureData : 'ctime'})#always use the most recent
    structures = q.all()
    
    if structures:
        return structures[-1][0]            
    else:
        return None

    
def fleur_calc_get_structure(calc_node):
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
    calc_uuids = []
    for calc in calcs:
        #print(calc)
        calc_uuids.append(calc.get_outputs_dict()['output_scf_wc_para'].get_dict()['last_calc_uuid'])
        #calc_uuids.append(calc['output_scf_wc_para'].get_dict()['last_calc_uuid'])
    #print(calc_uuids)
    
    all_corelevels = {}
    fermi_energies = {}
    bandgaps = {}
    all_atomtypes = {}  
    # more structures way: divide into this calc and reference calcs.
    # currently the order in calcs is given, but this might change if you submit
    # check if calculation pks belong to successful fleur calculations
    for uuid in calc_uuids:
        calc = load_node(uuid)
        if (not isinstance(calc, FleurCalc)):
            #raise ValueError("Calculation with pk {} must be a FleurCalculation".format(pk))
            # log and continue
            continue
        if calc.get_state() != 'FINISHED':
            # log and continue
            continue
            #raise ValueError("Calculation with pk {} must be in state FINISHED".format(pk))
        
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
        
        # TODO: maybe different, because it is prob know from before
        fleurinp = calc.inp.fleurinpdata
        structure = fleurinp.get_structuredata(fleurinp)            
        compound = structure.get_formula()
        #print compound
        fermi_energies[compound] = efermi
        bandgaps[compound] = bandgap
        all_atomtypes[compound] = atomtypes
        all_corelevels[compound] = corelevels
        #fermi_energies = efermi
        #bandgaps = bandgap
        #all_atomtypes = atomtypes
        #all_corelevels = corelevels
    
    return fermi_energies, bandgaps, all_atomtypes, all_corelevels
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
    from aiida.orm import Group
    from string import digits
    
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
            #print(message)
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
            #print(message)
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
    from aiida.orm import Group
    #from string import digits
    
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
            message = ('You have to provide a valid pk for a Group of' 
                      'parameters or a Group name. Reference key: "group".'
                      'given pk= {} is not a valid group'
                      '(or is your group name integer?)'.format(group_pk))
            #print(message)
            report.append(message)
    else:
        try:
            para_group = Group.get_from_string(group_name)
        except NotExistent:
            para_group = None
            message = ('You have to provide a valid pk for a Group of' 
                      'parameters or a Group name. Wf_para key: "para_group".'
                      'given group name= {} is not a valid group'
                      '(or is your group name integer?)'.format(group_name))
            #print(message)
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

'''
   def get_calcs_from_groups(self):
        """
        Extract the crystal structures and parameter data nodes from the given 
        groups and create calculation 'pairs' (stru, para).
        """
        wf_dict = self.inputs.wf_parameters.get_dict()
        #get all delta structure

        str_gr = wf_dict.get('struc_group', 'delta')
        
        try:
            group_pk = int(str_gr)
        except ValueError:
            group_pk = None
            group_name = str_gr
        
        if group_pk is not None:
            try:
                str_group = Group(dbgroup=group_pk)
            except NotExistent:
                str_group = None
                message = ('You have to provide a valid pk for a Group of' 
                          'structures or a Group name. Wf_para key: "struc_group".'
                          'given pk= {} is not a valid group'
                          '(or is your group name integer?)'.format(group_pk))
                #print(message)
                self.report(message)
                self.abort_nowait('I abort, because I have no structures to calculate ...')
        else:
            try:
                str_group = Group.get_from_string(group_name)
            except NotExistent:
                str_group = None
                message = ('You have to provide a valid pk for a Group of' 
                          'structures or a Group name. Wf_para key: "struc_group".'
                          'given group name= {} is not a valid group'
                          '(or is your group name integer?)'.format(group_name))
                #print(message)
                self.report(message)
                self.abort_nowait('I abort, because I have no structures to calculate ...')
                                
                
        #get all delta parameters
        para_gr = wf_dict.get('para_group', 'delta')
        
        if not para_gr:
            #waring use defauls
            message = 'COMMENT: I did recieve "para_group=None" as input. I will use inpgen defaults'
            self.report(message)
        
        try:
            group_pk = int(para_gr )
        except ValueError:
            group_pk = None
            group_name = para_gr 
        
        if group_pk is not None:
            try:
                para_group = Group(dbgroup=group_pk)
            except NotExistent:
                para_group = None
                message = ('ERROR: You have to provide a valid pk for a Group of' 
                          'parameters or a Group name (or use None for inpgen defaults). Wf_para key: "para_group".'
                          'given pk= {} is not a valid group'
                          '(or is your group name integer?)'.format(group_pk))
                #print(message)
                self.report(message)
                self.abort_nowait('ERROR: I abort, because I have no paremeters to calculate and '
                                  'I guess you did not want to use the inpgen default...')
        else:
            try:
                para_group = Group.get_from_string(group_name)
            except NotExistent:
                para_group = None
                message = ('ERROR: You have to provide a valid pk for a Group of' 
                          'parameters or a Group name (or use None for inpgen defaults). Wf_para key: "struc_group".'
                          'given group name= {} is not a valid group'
                          '(or is your group name integer?)'.format(group_name))
                #print(message)
                self.report(message)
                self.abort_nowait('ERROR: I abort, because I have no paremeters to calculate and '
                                  'I guess you did not want to use the inpgen default...')        

        # creating calculation pairs (structure, parameters)

        para_nodes = para_group.nodes
        n_para = len(para_nodes)
        stru_nodes = str_group.nodes
        n_stru = len(stru_nodes)
        if n_para != n_stru:
            message = ('COMMENT: You did not provide the same number of parameter'
                       'nodes as structure nodes. Is this wanted?')
            self.report(message)
        
        calcs = []
        for struc in stru_nodes:
            para = get_paranode(struc, para_nodes)
            calcs.append((struc, para))
        print calcs[:10]
        self.ctx.calcs_to_run = calcs
'''
