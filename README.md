# SBPEM_A
A Spherical Barotropic Primitive Equation Model on A grid

This Model was written by Wang Bin, and modified by Zhou Lilong at July 2018.

You can build up this model by a simple 'make' command, then a SBPEM.exe file will be generated.

The default compiler is gfortran, if you are going to use another one, modify configure.SBPEM as what you need.

Characteristics of this model:

Shallow water equations with IAP transformation;

Time integrate scheme: 
1. Consisitent Dissipation Operator(CDO in B.f90),
2. A new kind of 4th order Runge-Kutta(RK4.f90)
3. Predict-Correct(Predict-Correct.f90)

Antisymmetry Opterator(also the spatial discretization) in L.f90;

Split Scheme: 2nd order Conservative Split Pattern(CSP2.f90)

Parameter setting in module_para.f90, you can adjust the model resolution in this file, after changing the parameters,

you should recompile the model by following commad:

./clean

make

The forecasting result will be written into a netCDF file named 'output.nc'.

Zhou Lilong

National Meteorological Center of CMA

College of Earth and Planetary Sciences, UCAS

# Update log
V2.0
1. Add split scheme CSP2(2nd Order Conservative Split Pattern)
2. Add 4th order Improved Runge-Kutta and Predict-Correct as integrate scheme
