using Plots
using PyPlot
using LinearAlgebra
"""
 ~gergaud/ENS/Control/ODE/Matlab/test_gauss.m

 Auteurs:  Joseph GERGAUD
 Date:     nov. 2005
 Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
           2, rue Camichel 31071 Toulouse FRANCE
 Email:    gergaudenseeiht.fr
***************************************************************************

 Int�gration de l'equation differentiel de l'equation de Van der Pol
 ref: Hairer
"""

include("ode_gauss_v2.jl")
include("ode_gauss_v3.jl")
include("phi_vdp.jl")
include("d_phi_vdp.jl")
include("plot_sol.jl")
#rm("testgauss.txt")
#diary("testgauss.txt")
y0 = [2.008619860874843136; 0]
t0 = 0
tf = 6.663286859323130189
N0 = [120:60:1080; 1200:600:10800] # multiple de 12 pour avoir des nombres entier si on divise par 4 ou 3
N = N0
N = 10
#
# Gauss point fixe
N = 10
options = zeros(3)
options[1] = N
options[2] = 15
options[3] = 1e-6
T,Y,nphie,ifail,KK = ode_gauss_v2(phi_vdp,[t0 tf],y0,options)
#
println("Point fixe")
println("----------")
#println("[T Y]")
#println([T Y])
#println(nphie)
#println(ifail)
I = findall(ifail==-1)
ifail[I] .= options[2]

pyplot()
plt = Plots.plot(layout=(1,3))
plot_sol(plt,T,Y,"red","Gauss v2,N="*string(N))

# KK(:,cumsum(ifail+1))
#
# Gauss Newton
T,Y,nphie,ndphie,ifail,KK = ode_gauss_v3(phi_vdp,d_phi_vdp,[t0 tf],y0,options)
println("Newton")
println("------")
#println("[T Y]")
#println([T Y])
#println(nphie)
#println(ndphie)
#println(ifail)
#println(KK)
#

plot_sol(plt,T,Y,"green","Gauss v3,N="*string(N))
N = 200
options[1] = N

T,Y,nphie,ifail = ode_gauss_v2(phi_vdp,[t0 tf],y0,options)
plot_sol(plt,T,Y,"blue","Gauss v2,N="*string(N))