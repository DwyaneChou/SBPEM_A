include ./configure.SBPEM

OBJS = module_para.o      \
       module_array.o     \
       mesh.o             \
       haurwitz.o         \
       B.o                \
       CDO.o              \
       RK4.o              \
       Predict_Correct.o  \
       L.o                \
       fast_operator.o    \
       slow_operator.o    \
       CSP2.o             \
       integrator.o       \
       spatial_discrete.o \
       output_netCDF.o

all: EXE

EXE:  $(OBJS) main.o
	$(F90) -o $(EXENAME) $(OBJS) main.o \
	-L$(NETCDF)/lib -lnetcdf -lnetcdff $(OPT) $(FCFLAGS) $(CFLAGS)

module_para.o :
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) module_para.f90

module_array.o : module_para.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) module_array.f90

haurwitz.o :
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) haurwitz.f90

mesh.o : module_para.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) mesh.f90

L.o : module_para.o mesh.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) L.f90
	
fast_operator.o : module_para.o mesh.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) fast_operator.f90
	
slow_operator.o : module_para.o mesh.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) slow_operator.f90
	
B.o :
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) B.f90
	
CDO.o : module_para.o module_array.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) CDO.f90
	
RK4.o : module_para.o module_array.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) RK4.f90
	
Predict_Correct.o : module_para.o module_array.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) Predict_Correct.f90
  
integrator.o : module_para.o CDO.o RK4.o Predict_Correct.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) integrator.f90
  
spatial_discrete.o : module_para.o fast_operator.o slow_operator.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) spatial_discrete.f90
	
CSP2.o : module_para.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) CSP2.f90

output_netCDF.o : module_para.o mesh.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) output_netCDF.f90

main.o : $(OBJS)
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) main.f90


clean:
	rm *.o *.mod *.exe *.a
