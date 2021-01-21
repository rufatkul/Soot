# MOMIC MakeFile

module_name = momic_module
fort_files = momic_parameters.f90 momic.f90
F90FLAGS = F90FLAGS = -O0 -fPIC -m64 -g -traceback -fpp -fpe0 -ftrapuv -check bounds -check pointer -check all -DSERIAL


# Make File
$(module_name): $(fort_files)
	f2py -c $(fort_files) -m $(module_name) #--opt='-fbacktrace -fbounds-check'
