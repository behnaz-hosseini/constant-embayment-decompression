# constant-embayment-decompression
A MATLAB code (written in v. R2020a) for modeling diffusion of H2O and CO2 through rhyolitic melt embayments during constant, continuous magma
decompression. The example is set up to fit H2O and CO2 concentration gradients measured along a 90-um long quartz-hosted rhyolitic melt embayment that was
experimentally decompressed from 150 to 30 MPa at 780 C and a rate of 0.008 MPa/s. Please cite Hosseini et al. (2023), in revision at Geochemistry,
Geophysics, Geosystems.

Instructions for running code:

Download both emb_diffusion.m (master script) and diffusion_function.m (function) into a folder.

emb_diffusion.m is the main code where all input parameters are defined and the misfit is calculated. On Lines 20-26, specify input parameters, including
initial pressure (MPa), embayment length (um), initial H2O (wt. %), melt density (kg/m^3), temperature (C), and uncertainties. On lines 33-35: enter H2O
+/- CO2 concentrations (wt. % and ppm, respectively), as well as distance along the measured profile (um). All arrays go from embayment outlet to the
interior. On lines 39-41: specify the range of decompression rates (MPa/s), fragmentation pressures (MPa), and bubble radii to iterate through.

Diffusion_function_1D_constant.m is the function that performs the diffusion calculation. User does not need to modify.
