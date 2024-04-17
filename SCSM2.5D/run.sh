gfortran -mcmodel=medium  -c C_DF.f90 Grid_Model.f90 Interp.f90 MATRIX_YYXZ.f90 MS_DF.f90 Gauss_Quad.f90 Viscoelastic.f90 Viscoelastic2.f90 
gfortran -mcmodel=medium  C_DF.o Grid_Model.o Interp.o MATRIX_YYXZ.o MS_DF.o Gauss_Quad.o Viscoelastic.o Viscoelastic2.o MainOMP.f90 -fopenmp 
export omp_set_num_threads=2 

./a.out

