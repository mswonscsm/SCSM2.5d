# SCSM2.5d
SCSM2.5D is a 2.5D seismic wave modeling program that uses a subdomain Chebyshev spectral finite difference method. 
The 2.5D modeling technique produces a 3D wavefield using a 2D geological model, making it ideal for seismic line surveys. 
The program can handle various modeling scenarios including viscoacoustic, isotropic viscoelastic, anisotropic viscoelastic (VTI, ORT, and TTI), solid and water free surfaces, and water-solid interfaces.

# Installment
SCSM2.5d OpenMP version: Compile MainOMP.f90, C_DF.f90, Gauss_Quad.f90, Grid_Model.f90, interp.f90, MATRIX_YYXZ.f90, MS_DF.f90, Viscoelastic.f90, and Viscoelastic2.f90. 

SCSM2.5d MPI/OpenMP version: Compile MainMPIOMP.f90, C_DF.f90, Gauss_Quad.f90, Grid_Model.f90, interp.f90, MATRIX_YYXZ.f90, MS_DF.f90, Viscoelastic.f90, and Viscoelastic2.f90.
By substituting MainOMP.f90 with MainMPIOMP.f90, the program can be converted to an MPI/OpenMP version that supports fully-parallel computation, utilizing cores equal to the number of wavenumber samples. 
This enhancement reduces computation time to less than 1.5 times that of 2D modeling. 

Windows: If you don't have Fortran compiler, you try the 'target.exe' with different input files (2.5Dseis_SCSM.inp and relaxation_time.inp) in Example folder.

Linux: run 'run.sh'.

HPC: run 'runOMP.sh' or 'runMPIOMP.sh'.
