#
# Intégration de L'équation différentielle considérée est l'équation de Van der Pol
# ref: Hairer
#

using Plots
using PyPlot
using LinearAlgebra

include("solvers/ode_gauss_newton.jl")
include("solvers/ode_gauss_pf.jl")
include("solvers/ode_gauss.jl")
include("solvers/ode_euler.jl")
include("solvers/ode_runge.jl")
include("solvers/ode_rk41.jl")
include("solvers/ode_rk42.jl")
include("solvers/ode_heun.jl")
include("plot_sol.jl")
include("functions/d_fun_vdp.jl")
include("functions/fun_vdp.jl")

y0 = [2.008619860874843136; 0]
t0 = 0
tf = 6.663286859323130189
N = 25
closeall()
pyplot() # utiliser Pyplot comme backend

# Gauss point fixe
options = zeros(3)
options[1] = N
options[2] = 15
options[3] = 1e-12
T,Y,nphie,ifail = ode_gauss_pf(fun_vdp,[t0 tf],y0,options)
#
println("Point fixe")
println("----------")
println("[T Y]")
display("text/plain",[T Y])
println("nphie")
display("text/plain",nphie)
println("ifail")
display("text/plain",ifail)
#I = findall(ifail==-1)
#ifail[I] .= options[2]
plt = Plots.plot(layout=(3))
plot_sol(plt,T,Y,"red","Gauss v2,N="*string(N))

# Gauss Newton
T,Y,nphie,ndphie,ifail = ode_gauss_newton(fun_vdp,d_fun_vdp,[t0 tf],y0,options)
println("Newton")
println("------")
println("[T Y]")
display("text/plain",[T Y])
println("nphie")
display("text/plain",nphie)
println("ndphie")
display("text/plain",ndphie)
println("ifail")
display("text/plain",ifail)
plot_sol(plt,T,Y,"darkorange","Newton,N="*string(N))

T,Y = ode_euler(fun_vdp,[t0 tf],y0,N)
println("Euler")
println("------")
println("[T Y]")
display("text/plain",[T Y])
plot_sol(plt,T,Y,"blue","Euler,N="*string(N))

T,Y = ode_runge(fun_vdp,[t0 tf],y0,N)
println("Runge")
println("------")
println("[T Y]")
display("text/plain",[T Y])
plot_sol(plt,T,Y,"black","Runge,N="*string(N))

T,Y = ode_heun(fun_vdp,[t0 tf],y0,N)
println("Heun")
println("------")
println("[T Y]")
display("text/plain",[T Y])
plot_sol(plt,T,Y,"magenta","Heun,N="*string(N))

T,Y = ode_rk41(fun_vdp,[t0 tf],y0,N)
println("RK41")
println("------")
println("[T Y]")
display("text/plain",[T Y])
plot_sol(plt,T,Y,"green","RK41,N="*string(N))

T,Y = ode_rk42(fun_vdp,[t0 tf],y0,N)
println("RK42")
println("------")
println("[T Y]")
display("text/plain",[T Y])
plot_sol(plt,T,Y,"cyan","RK42,N="*string(N))

#N = 200
#options[1] = N
#
#T,Y,nphie,ifail = ode_gauss_pf(fun_vdp,[t0 tf],y0,options)
#plot_sol(plt,T,Y,"blue","Gauss v2,N="*string(N))