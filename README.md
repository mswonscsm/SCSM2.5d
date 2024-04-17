# SCSM2.5d
SCSM2.5D is a 2.5D seismic wave modeling program that uses a subdomain Chebyshev spectral finite difference method. 
The 2.5D modeling technique produces a 3D wavefield using a 2D geological model, making it ideal for seismic line surveys. 
The program can handle various modeling scenarios including viscoacoustic, isotropic viscoelastic, anisotropic viscoelastic (VTI, ORT, and TTI), solid and water free surfaces, and water-solid interfaces.

# Installment
* SCSM2.5d OpenMP version: Compile __MainOMP.f90, C_DF.f90, Gauss_Quad.f90, Grid_Model.f90, interp.f90, MATRIX_YYXZ.f90, MS_DF.f90, Viscoelastic.f90, Viscoelastic2.f90__. 
* SCSM2.5d MPI/OpenMP version: Compile __MainMPIOMP.f90, C_DF.f90, Gauss_Quad.f90, Grid_Model.f90, interp.f90, MATRIX_YYXZ.f90, MS_DF.f90, Viscoelastic.f90, and Viscoelastic2.f90__.

By substituting MainOMP.f90 with MainMPIOMP.f90, the program can be converted to an MPI/OpenMP version that supports fully-parallel computation, utilizing cores equal to the number of wavenumber samples. 
This enhancement reduces computation time to less than 1.5 times that of 2D modeling. 

* Windows: If you don't have Fortran compiler, you can try the __'target.exe'__ with different input files __(2.5Dseis_SCSM.inp and relaxation_time.inp)__ in Example folder.
However, you need to optimize absorbing layer, 2.5-D setting, core number for parallel computing, etc in Fortran code for your case. So, just try __target.exe__ for program test only. 
* Linux: run __'run.sh'__.
* HPC: run __'runOMP.sh' or 'runMPIOMP.sh'__.

# Input data
SCSM2.5d has two input data, and details are explained in manual, Doc.
* __2.5Dseis_SCSM.inp__
* __relaxation_time.inp__

# Output data
* Waveform figure: __fort.xxx, Xgrid.out, Zgrid.out__ 
You can generate the figure in __Matlab Figure Waveform__.
__fort.xxx__ is generated only in 2D modeling.

* Seismogram: __rec_real_x.out__
You can make seismogram in __Matlab Figure Seismogram__ and __Matlab Figure SeismicLine__ (for massive and complex seismograms).

# Contact
* Moosoo Won, PhD in Earth Sciences, Khalifa University of Science and Technology, 100058280@ku.ac.ae, merccer999@gmail.com
* Dr.Bing Zhou, Associate Professor, Earth Sciences, Khalifa University of Science and Technology, bing.zhou@ku.ac.ae
