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

include("ode_gauss_v3.jl")
include("ode_gauss_pf.jl")
include("ode_gauss.jl")
include("ode_euler.jl")
include("ode_runge.jl")
include("ode_rk41.jl")
include("ode_rk42.jl")
include("ode_heun.jl")
include("plot_sol.jl")
include("d_fun_vdp.jl")
include("fun_vdp.jl")
#rm("testgauss.txt")
#diary("testgauss.txt")
y0 = [2.008619860874843136; 0]
t0 = 0
tf = 6.663286859323130189
N = 25
#
# Gauss point fixe
options = zeros(3)
options[1] = N
options[2] = 15
options[3] = 1e-12
T,Y,nphie,ifail,KK = ode_gauss_pf(fun_vdp,[t0 tf],y0,options)
#
println("Point fixe")
println("----------")
#println("[T Y]")
#println([T Y])
#println(nphie)
#println(ifail)
#I = findall(ifail==-1)
#ifail[I] .= options[2]
closeall()
pyplot()
plt = Plots.plot(layout=(1,3))
plot_sol(plt,T,Y,"red","Gauss v2,N="*string(N))

# KK(:,cumsum(ifail+1))
#
# Gauss Newton
T,Y,nphie,ndphie,ifail,KK = ode_gauss_v3(fun_vdp,d_fun_vdp,[t0 tf],y0,options)
println("Newton")
println("------")

plot_sol(plt,T,Y,"darkorange","Newton,N="*string(N))

T,Y = ode_euler(fun_vdp,[t0 tf],y0,N)
println("Euler")
println("------")
plot_sol(plt,T,Y,"blue","Euler,N="*string(N))

T,Y = ode_runge(fun_vdp,[t0 tf],y0,N)
println("Runge")
println("------")

plot_sol(plt,T,Y,"black","Runge,N="*string(N))

T,Y = ode_heun(fun_vdp,[t0 tf],y0,N)
println("Heun")
println("------")
plot_sol(plt,T,Y,"magenta","Heun,N="*string(N))

T,Y = ode_rk41(fun_vdp,[t0 tf],y0,N)
println("RK41")
println("------")
plot_sol(plt,T,Y,"green","RK41,N="*string(N))

T,Y = ode_rk42(fun_vdp,[t0 tf],y0,N)
println("RK42")
println("------")
plot_sol(plt,T,Y,"cyan","RK42,N="*string(N))

#N = 200
#options[1] = N
#
#T,Y,nphie,ifail = ode_gauss_pf(fun_vdp,[t0 tf],y0,options)
#plot_sol(plt,T,Y,"blue","Gauss v2,N="*string(N))