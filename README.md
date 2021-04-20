# gFDM
Ghost-filling Finite Difference Method

gFDM Implementation
gFDM implementation is divided in 4 high level functions and 32 low level private functions. All functions start with the “gfdm_” prefix text. Figure 1 illustrates the pipeline for the high-level functions, starting from the NET preprocessing pipeline. The method supports the inclusion of CTI maps that can handle isotropic and anisotropic conductivity definitions. The pipe line is divided in 3 steps, namely the head model “gfdm_prepare_headmodel” routine to calculate the stiffness matrix, the “gfdm_prepare_elecs” routine to check the electrodes positions and setup the lead-pair calculations, finally, the forward solutions routines “gfdm_precalculate_leads” and “gfdm_calculate_pots” calculates the reciprocity lead pair potentials and the output potentials for a given source space respectively.


High level function descriptions

 ![image](https://user-images.githubusercontent.com/49439997/115318697-e1e9aa80-a143-11eb-9e6d-439fe6368606.png)

Figure 1: Net headmodelling - gFDM pipeline implementation

gfdm_prepare_headmodel calculates the sparse stiffness matrix to solve the linear system. The function also finds the boundary box enclosing the head volume, labeling the voxel positions excluding the air. The output is a structure including the stiffness matrix and the boundary box array for the head volume. The gFDM method allows arbitraries voxel-size, a rectangular 1mm^3 voxel size is not mandatory, and the user can define not rectangular voxels. Also, using a down-sampling function (as “ft_volumedownsample”) the forward calculations can be drastically reduced using a low-resolution head volume. In the cfg structure, one can choose between a labeled map, where the cond_value array corresponds to the isotropic conductivities for the N layers in the MRI segmentation, or a CTI including a conductivity tensor definition for every single voxel in the head volume.

gfdm_prepare_headmodel code parameters:

INPUT: ( cfg ) The cfg argument is a structure containing:
cfg.downsample               = downsample value (integer);
cfg.resolution               = voxel size in mm; In mm
cfg.segmentation             = structural segmented MRI;
cfg.conductivity             = structure holding the conductivity paremeters;
    conductivity.tissuelabel = array holding the tissue names;
    conductivity.cond_value  = vector holding the tissue conductivity values;
cfg.conductivity_map         = (optional) structural CTI isotropic map;
cfg.gm_idx                   = label value of the grey matter;
vol  = gfdm_prepare_headmodel(cfg); 

OUTPUT: ( vol ) structure containing the stiff matrix
vol.box          = structure holding the boundary box;
vol.stiff        = sparse stiffness matrix
vol.c_idx        = 3D array holding the non-zero index.
vol.Cidx         = 1D array holding the non-zero index.
vol.Nz           = Number of non-zeros.

EXAMPLE: example using a CTI map
cfg                  = [];
cfg.resolution       = [1 1 1];
cfg.segmentation     = segmented_mri;
cfg.conductivity_map = CTI;
cfg.gm_idx           = 2;
vol = gfdm_prepare_headmodel(cfg);


gfdm_prepare_elecs ensures that the already aligned electrodes are placed over the scalp surface, projecting the electrode positions across a line with the center in the middle point of the GM. The electrodes must be placed over scalp voxel positions, otherwise, the lead-pair source will not have a valid position in the right-hand side vector. 

gfdm_prepare_elecs code parameters:

INPUT: ( cfg ) The cfg argument is a structure containing:
cfg.mri           a fieldtrip MRI structure holding an alrready aligned structural     MRI.
cfg.electrodes    a fieldtrip electrode structure containing the aligned electrodes.
cfg.vol           head model structure.

OUTPUT: ( elec ) structure containing the FDM electrodes
  elec_fdm.lead    array holding the lead-pair distribution
  elec_fdm.leadp   vector holding the lead-pair distribution
  elec_fdm.cent    midle point of the segmented GM
  elec_fdm.OK      boolean flag indicating 

EXAMPLE:
cfg               = [];
cfg.mri           = mri_aligned;
cfg.electrodes    = electrodes_aligned;
cfg.vol           = vol;
elec              = gfdm_prepare_elecs( cfg );

gfdm_precalculate_leads calculates lead-pair potentials for a given set of electrode positions. The calculations are performed using the reciprocity theorem for finite differences. This routine takes about 5 minutes per lead-pair (par of electrodes), using the number of electrodes minus one to complete the task. Each lead-pair is solved as a single linear system using iLU factorization preconditioner and BiCG-Stabilized solver. Time can be larger than 5 minutes if the stiff matrix contains anisotropic conductivity values. Eg: for a set of 70 electrodes, 1mm^3 voxel size, realistic head model, with +5million Non-zero voxels, this routine can yield +6h of calculations. in addition, the amount of memory needed to allocate the output is directly proportional to the amount of non-zero voxels. For the example mentioned the amount of memory used for the output is +1.4Gb. Although it is possible to calculate the potentials just for the GM area, or even for a given source distribution, drastically reducing the amount of memory needed to allocate the output.

gfdm_precalculate_leads code parameters:

INPUT: ( cfg ) The cfg argument is a structure containing:
cfg.vol            head model structure, output -> gfdm_prepare_headmodel
cfg.electrodes     electrodes structure, output -> gfdm_prepare_elecs
cfg.mri            structural aligned MRI
cfg.solver         a structure holding the linear solver paremeters
    solver.tol     tolerance error;
    solver.maxit   maximun number of iterations
cfg.gmask          a string that allows the reduction of the stored potentials 
                   in the head volume. This parameter can be  'source', 'gm', or 'all'.
                   'all' -> will return the potentials for the entire volume.
                   'gm' -> will return the potentials just for the GM 
                   tissue volume, allowing not only a high reduction of 
                   the memory needed for the precalculated leal-pairs but 
                   also to select any position/orientation for sources placed. inside the GM.
                   'source' -> will return the potentials needed for a grid 
                   source distribution obtained using ft_prepare_sourcemodel. 
                   Eg: grid = ft_prepare_sourcemodel(cfg)
cfg.grid           a structur containing a regular 3D source grid
                   obtained using ft_prepare_sourcemodel. Just for
                   gmask == 'source'.
                   Eg: grid = ft_prepare_sourcemodel(cfg)

OUTPUT: ( lead ) structure containing the precalculated lead-pairs
lead.vol           head model structure, output -> gfdm_prepare_headmodel
lead.electrodes    electrodes structure, output -> gfdm_prepare_elecs  
lead.mri           structural aligned MRI
lead.c_idx         3D array holding the Non-zero indexation
lead.solver        structur holding the performance parameters for the
                   linear solver in every lead-pair
     solver.tleads    vector containing the spend time for every lead-pair
     solver.fla       vector containing the final flags for the BiCG-Stabilized solver
     solver.rra       vector holding the final relative residual errors
     solver.ita       vector holding the number of iterations needed to solve every   lead-pair.
     solver.rva       cell array containing the relative error progression for every lead-pair.
lead.indexes          sparse vector containing the Non-zero potentials indexes
lead.pots             sparce matrix (number of voxels x 3) holding the x,
                      y and z potentials in the head volume.

EXAMPLE: example using already corregistered electrodes and aligned MRI data
cfg               = [];
cfg.vol           = vol;
cfg.elec          = elec_FDM;
cfg.mri           = mri_aligned;
cfg.solver.tol    = 1e-10;
cfg.solver.maxit  = 700;
cfg.gmask = 'gm';
lead = gfdm_precalculate_leads( cfg );


gfdm_calculate_pots calculates the leadfield matrix for a given grid-source space and precalculated lead-pairs potentials (gfdm_precalculate_leads output).

gfdm_calculate_pots code parameters:

INPUT: ( lead, cfg ). lead is the ouput of the gfdm_precalculate_leads function. The cfg argument is the same structure for the gfdm_precalculate_leads funtion imput containing:
cfg.vol            head model structure, output -> gfdm_prepare_headmodel;
cfg.electrodes     electrodes structure, output -> gfdm_prepare_elecs;
cfg.mri            structural aligned MRI;
cfg.solver         structure holding the linear solver paremeters;
cfg.gmask          string that allows the reduction of the stored potentials;
cfg.grid           structur containing source position;

OUTPUT: ( leadfield ) structure holding the leadfield potentials for a given 
grid for source positions.
leadfield.inside       logical vector holding flags for the source positions
                       inside the grey matter;
leadfield.pos          2D array holding the source positions;
leadfield.leadfield    cell array holding the leadfiled potentials

EXAMPLE
lead -> output -> gfdm_precalculate_leads
cfg  -> input  -> gfdm_precalculate_leads
leadfield = gfdm_calculate_pots(lead, cfg);

