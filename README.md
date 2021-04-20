# gFDM
Ghost-filling Finite Difference Method


1.	gFDM Implementation
gFDM implementation is divided in 4 high level functions and 32 low level private functions. All functions start with the “gfdm_” prefix text. Figure 1 illustrates the pipeline for the high-level functions, starting from the NET preprocessing pipeline. The method supports the inclusion of CTI maps that can handle isotropic and anisotropic conductivity definitions. The pipe line is divided in 3 steps, namely the head model “gfdm_prepare_headmodel” routine to calculate the stiffness matrix, the “gfdm_prepare_elecs” routine to check the electrodes positions and setup the lead-pair calculations, finally, the forward solutions routines “gfdm_precalculate_leads” and “gfdm_calculate_pots” calculates the reciprocity lead pair potentials and the output potentials for a given source space respectively.


1.1.	High level function descriptions

 ![image](https://user-images.githubusercontent.com/49439997/115318697-e1e9aa80-a143-11eb-9e6d-439fe6368606.png)

Figure 1: Net headmodelling - gFDM pipeline implementation

gfdm_prepare_headmodel calculates the sparse stiffness matrix to solve the linear system. The function also finds the boundary box enclosing the head volume, labeling the voxel positions excluding the air. The output is a structure including the stiffness matrix and the boundary box array for the head volume. The gFDM method allows arbitraries voxel-size, a rectangular 1mm^3 voxel size is not mandatory, and the user can define not rectangular voxels. Also, using a down-sampling function (as “ft_volumedownsample”) the forward calculations can be drastically reduced using a low-resolution head volume. In the cfg structure, one can choose between a labeled map, where the cond_value array corresponds to the isotropic conductivities for the N layers in the MRI segmentation, or a CTI including a conductivity tensor definition for every single voxel in the head volume.
