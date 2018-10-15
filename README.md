# MLEM-source-reconstruction

A code to reconstruct linax X-ray sources based on crossplane and inplane profile measurements in air. 

The srcrec_main is the main script to run. Other scripts needed with it: 
srcrec_errors: reconstructes the source multiple times and provides an estimate of the total ucnertainty,
RecSource: Reconstructs the source based on input field parameters and system matrix,
mlem: the mlem iterative algorithm as applied in one full reconstruction,
ExtrSystemMat : extracts the system matrix for the given linac gometry

INPUT variables: x (mm): off-axis positions of dose measurements (-8 mm -> 8mm, step 0.2 mm), dose: relative dose profile,
dose_std: st.dev of relative dose profile measurements, pro: 'cro' or 'in' for crossplane or inplane,
orientation, N: number of reconstructions to run for total uncertainty estimation (typically N = 30 or higher)

DATA FILES
PSF.mat: The PSF kernel needed for the profile deconvolution step (this needs to be in the same directory as the script files)
example_sources.mat: 3 example sets of crossplane and inplane profile measurements. Can be used as input in srcerec_main

OUTPUT: A: System matrix, field_opt (mm) : optimum field size selected, 
snrec: reconstructed source, snrec_std: st dev of rec source,
FWHMrec (mm) : FWHM of reconstructed source, TWHMrec (mm) : TWHM of rec source, 
FWHMstd (mm) : std of FWHM of reconstructed source, TWHMstd (mm) : std of TWHM of rec source, 
pro_err90_10 : mean absolute error between input and rec profile in the 90-10 % region, 
src_RMSE100_10 : rmse between the rec source and a gaussian fit in the 100 - 10 % region

The code outputs also a figure exhibiting (a) the reconstructed source and a gaussian fit, 
(b) measured and reconstructed profiles and (c) the field opimization cost function. 

A detailed description of the method and validation of the code can be found in this article (open access): 

http://iopscience.iop.org/article/10.1088/0031-9155/61/3/1078/meta
