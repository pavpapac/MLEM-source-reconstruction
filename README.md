# MLEM-source-reconstruction

A code to reconstruct linac X-ray sources based on photon fluence profile measurements in air. 

MAIN code
    
    srcrec_main

    INPUT    
    x (mm) : off-axis positions of dose measurements (see example_sources)
    dose : a 1-d array of a relative dose profile measurements (see example_sources)
    dose_std : a 1-d array of st.dev of relative dose profile measurements (see example_sources)
    
    OUTPUT 
    A: The linac geometry system matrix,
    field_opt (mm) : the optimum field size selected, 
    snrec: the reconstructed source fluence profile, 
    snrec_std: the st dev of reconstructed source fluence profile,
    FWHMrec (mm) : FWHM of reconstructed source,
    FWHMstd (mm) : std of FWHM of reconstructed source,
    TWHMrec (mm) : TWHM of rec source, 
    TWHMstd (mm) : std of TWHM of rec source, 
    pro_err90_10 : mean absolute error between input and reconstructed profile in the 90-10 % dose region, 
    src_RMSE100_10 : rmse between the rec source and a gaussian fit in the 100-10 % fluence region
    
The code outputs also a figure exhibiting (a) the reconstructed source and a gaussian fit, 
(b) the measured and reconstructed profile and (c) the field opimization cost function.

SCRIPTS

    srcrec_errors: reconstructs the source multiple times and provides an estimate of the total ucnertainty,
    RecSource: Reconstructs the source based on input field parameters and system matrix,
    MLEM: the mlem iterative algorithm as applied in one full iteration,
    ExtrSystemMat : extracts the system matrix for the given linac gometry
    
DATA FILES

    PSF.mat: The PSF kernel needed for the profile deconvolution step. This needs to be in the same directory as the srcrec_main.
    example_sources.mat: 3 example sets of crossplane and inplane fluence profiles and st. dev. 


The source reconstruction process (summary)

The experimental measurements consist of direct measurements of the photon fluence profile in air using radiochromic films.Films are positioned beneath a 2 mm Pb foil build-up material. The Pb foil assists in reducing the profile blurring due to the non-zero electron range and photon scattering. In order to eliminate any residual blurring, a MC pre-calculated deconvolution kernel (PSF) is directly
applied on the dose profiles.  In order to minimize backscattering events a 5 cm styrofoam block is placed beneath the films. The
measurement plane is set at a Source to Surface distance (SSD) of 105 cm. The field size is set to 0.5 x 0.5 cm2 as defined by the secondary collimators, in order to reduce the impact of scatter sources and other parameters on the dose profile.
After the photon fluence profile has been extracted, it can be used as a direct input to the iterative reconstruction algorithm, based on the maximum likelihood expectation maximization (MLEM) algorithm. Firstly, the linac geometry is modeled
by the system matrix which maps each source pixel to the measurement plane. The procedure starts with an arbitrary assumption of the source distribution, typically a uniform distribution. Photons are then projected on the measurement plane by multiplying the system matrix with the source distribution array. Fluence profile correction per pixel are then extracted as the ratio of measured to expected photon  fluence. The corrections are re-normalized and ray-traced backwards to the source plane. The source distribution shape is updated by applying the corrections pixel-by-pixel. The algorithm continues iteratively until the FWHM of the source has not varied
more than 2 % in the last 30 iterations. A caveat on the previous procedure is that the actual collimator field size at the time of measurement is not known, For such narrow field sizes, even a small variation of the jaw positioning will significantly impact the source occlusion effect and thus the system matrix. In order to address this issue the jaw positions varied from 0.5 to 0.6 cm with a step 0.01 cm and the reconstruction was repeated for each value. The jaw position that resulted in the lowest mean local error between
calculated and measured dose profiles was selected as the optimum choice.


A detailed description of the method and validation of the code can be found in this article (open access): 

http://iopscience.iop.org/article/10.1088/0031-9155/61/3/1078/meta
