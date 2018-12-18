Fleur core-hole workflow
------------------------

Class name, import from:
  ::

    from aiida_fleur.workflows.corehole import fleur_corehole_wc
    #or 
    WorkflowFactory('fleur.corehole')

Description/Purpose
^^^^^^^^^^^^^^^^^^^

  Workflow for core-hole calculations, for the purpose of calculation electron binding energies.
  
  Runs a supercell calculation on the original structure and compares its total energy to 
  systems with a single core-hole. The atom with the corehole will be centered, and the symmetry of the cell will be broken.
  Specification of the core-hole are given in the input parameters.
  u.a. for what atoms, elements and corelevels corehole calculations should be performed.
  The corehole typ can be ``charge or valence``. in the charge case the electron in removed.
  In the valence case a core electron is transfered to the highest non full valence state as starting point.
  The workflow also allows for partial coreholes of any charge and spinpolarized calculations are the default.
  
  Converges the charge density and/or the total energy of a given system, 
  or stops because the maximum allowed retries are reached.
    
  This workflow manages none to one inpgen calculation and two to a lot of Fleur calculations.
  It is an quite advanced workflow which can spawn a lot of large scf workchains as subworkflows.
  
  .. note::

    
Input nodes
^^^^^^^^^^^
  * ``fleur`` (*aiida.orm.Code*): Fleur code using the ``fleur.fleur`` plugin
  * ``inpgen`` (*aiida.orm.Code*): Inpgen code using the ``fleur.inpgen`` plugin
  * ``wf_parameters`` (*ParameterData*, optional): Some settings of the workflow behavior (e.g. method, hole charge, atom, corelevel, ...)
  
  * ``structure`` (*StructureData*, path 1): Crystal structure data node.
  * ``calc_parameters`` (*str*, optional): specific FLAPW parameters to the corresponding structure
    
  * ``fleurinp`` (*FleurinpData*, path 2): Label of the workflow
  * ``remote_data`` (*RemoteData*, optional): The remote folder of the (converged) calculation whose output potential is used as input for the DOS run
  * ``options`` (*ParameterData*, optional): option 
# * ``settings`` (*ParameterData*, optional): special settings for Fleur calculations, will be given like it is through to calculationss.
    
Returns nodes
^^^^^^^^^^^^^
  * ``output_scf_wc_para`` (*ParameterData*): Information of workflow results like success, last result node, list with convergence behavior

  * ``fleurinp`` (*FleurinpData*) Input node used is retunred.
  * ``last_fleur_calc_output`` (*ParameterData*) Output node of last Fleur calculation is returned.

  
Default input
^^^^^^^^^^^^^

  * ``wf_parameters``: {
            'method' : 'valence', # what method to use, default for valence to highest open shell
            'hole_charge' : 1.0,       # what is the charge of the corehole? 0<1.0
            'atoms' : ['all'],           # coreholes on what atoms, positions or index for list, or element ['Be', (0.0, 0.5, 0.334), 3]
            'corelevel': ['all'],        # coreholes on which corelevels [ 'Be1s', 'W4f', 'Oall'...]
            'supercell_size' : [2,1,1], # size of the supercell [nx,ny,nz]
            'para_group' : None,       # use parameter nodes from a parameter group
            #'references' : 'calculate',# at some point aiida will have fast forwarding
            'relax' : False,          # relax the unit cell first?
            'relax_mode': 'Fleur',    # what releaxation do you want
            'relax_para' : 'default', # parameter dict for the relaxation
            'scf_para' : 'default',    # wf parameter dict for the scfs
            'same_para' : True,        # enforce the same atom parameter/cutoffs on the corehole calc and ref
            'resources' : {"num_machines": 1},# resources per job
            'max_wallclock_seconds' : 6*60*60,    # walltime per job
            'queue_name' : '',       # what queue to submit to
            'serial' : True,           # run fleur in serial, or parallel?
            #'job_limit' : 100          # enforce the workflow not to spawn more scfs wcs then this number(which is roughly the number of fleur jobs)
            'magnetic' : True          # jspins=2, makes a difference for coreholes
            }
Layout
^^^^^^
  .. figure:: /images/Workchain_charts_corehole_wc.png
    :width: 50 %
    :align: center

Database Node graph
^^^^^^^^^^^^^^^^^^^
  .. code-block:: python
    
    from aiida_fleur.tools.graph_fleur import draw_graph
    
    draw_graph(30528)
    
  .. figure:: /images/corehole_si_30528.pdf
    :width: 100 %
    :align: center
        
Plot_fleur visualization
^^^^^^^^^^^^^^^^^^^^^^^^
  Currently there is no visualization directly implemented for plot fleur.
  Through there are construct and plot spectra method from Binding energies and core level shifts in 
  ``masci-tools/vis/plot_methods.py``

Example usage
^^^^^^^^^^^^^
  .. include:: ../../../../examples/tutorial/workflows/tutorial_submit_corehole_wc.py
     :literal:

     
Output node example
^^^^^^^^^^^^^^^^^^^
  .. include:: /images/corehole_wc_outputnode.py
     :literal:

Error handling
^^^^^^^^^^^^^^
  Still has to be documented