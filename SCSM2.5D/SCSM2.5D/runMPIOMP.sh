#!/bin/bash
#SBATCH --nodes=3
#SBATCH --ntasks=78

#SBATCH --job-name=M_W_80
#SBATCH --time=48:00:00
#SBATCH --partition=prod
#SBATCH --account=kunf0069
#SBATCH --output=out.%j
#SBATCH --error=err.%j
 
module purge
module load gcc/9.4
module load openmpi/4.1


mpif90 -c C_DF.f90 Grid_Model.f90 Interp.f90 MATRIX_YYXZ.f90 MS_DF.f90 Gauss_Quad.f90 Viscoelastic.f90 Viscoelastic2.f90 
mpif90 C_DF.o Grid_Model.o Interp.o MATRIX_YYXZ.o MS_DF.o Gauss_Quad.o Viscoelastic.o Viscoelastic2.o MainMPIOMP.f90 -fopenmp

#ulimit -s unlimited
#export OMP_STACKSIZE=800g
export omp_set_num_threads=1 

mpirun -np 78 ./a.out | tee info10.log
