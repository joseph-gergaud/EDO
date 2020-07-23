#
# Intégration de l'équation differentiel de l'équation de Van der Pol
# ref: Hairer
#
using Plots

include("functions/fun_toto.jl")
include("solvers/ode_euler.jl")
include("plot_sol.jl")

y0 = [2.008619860874843136;  0]
t0 = 0
tf = 6.663286859323130189

T,Y = ode_euler(fun_toto,[t0 tf],y0,5)
plt = Plots.plot(layout=(1,3))
plot_sol(plt,T,Y,"blue","euler")
Plots.plot!(T*ones(1,2),Y)
