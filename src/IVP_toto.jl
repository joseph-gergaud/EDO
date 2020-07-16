using PyPlot

include("phi_toto.jl")
include("euler.jl")
#
# ~gergaud/ENS/Control/ODE/Matlab/IVP_toto.m
#
# Auteurs:  Joseph GERGAUD
# Date:     nov. 2005
# Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
#           2; rue Camichel 31071 Toulouse FRANCE
# Email:    gergaud@enseeiht.fr
#***************************************************************************
#
# Intégration de l"equation differentiel de l"equation de Van der Pol
# ref: Hairer
#
y0=[2.008619860874843136;  0]
t0=0
tf=6.663286859323130189
T,Y = ode_euler(phi_toto,[t0 tf],y0,5)
plot(T*ones(1,2),Y)
