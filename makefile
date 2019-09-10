#Generic makefile for fortran90


FC := ifort  

OBJ := dihedral_correlation.o\
       main.o

%.o : %.f90
	$(FC) -c $< 

all : exe

exe :  $(OBJ) 
	$(FC) -o $@  $^ 

main.o: dihedral_correlation.o

clean :
	rm *.o *.mod exe -r qm_dih_cor_dat -r qm_dih_d_value -r mm_dih_cor_dat -r mm_dih_d_value
