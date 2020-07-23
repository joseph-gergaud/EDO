"""
 Int�gration de l'equation differentiel de l'equation de Van der Pol
 ref: Hairer
"""

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

pause(text) = (print(stdout, text); read(stdin, 1); nothing)

y0 = [2.00861986087484313650940188;  0]
t0 = 0
tf = 6.6632868593231301896996820305
N0 = 10
N = N0
Ifail = zeros(length(N))
options = zeros(3)
#
# Euler

for i=1:length(N)
     println("Euler")
     T,Y = ode_euler(fun_vdp,[t0 tf],y0,N[i])
     display("text/plain",[T Y])
     display(Y)     
     pause("tapez entrée pour continuer")     
     plt = Plots.plot(layout=(1,2))
     Plots.plot!(T,Y[:,1],subplot=1,label="Euler")
     Plots.plot!(T,Y[:,2],subplot=2,label="Euler")

     println("Runge")
     T,Y=ode_runge(fun_vdp,[t0 tf],y0,N[i])
     display("text/plain",[T Y])
     Plots.plot!(T,Y[:,1],subplot=1,label="Runge")     
     Plots.plot!(T,Y[:,2],subplot=2,label="Runge")

     println("Heun")
     T,Y=ode_heun(fun_vdp,[t0 tf],y0,N[i])
     display("text/plain",[T Y])
     Plots.plot!(T,Y[:,1],subplot=1,label="Heun")     
     Plots.plot!(T,Y[:,2],subplot=2,label="Heun")

     println("RK41")
     T,Y=ode_rk41(fun_vdp,[t0 tf],y0,N[i])
     display("text/plain",[T Y])
     Plots.plot!(T,Y[:,1],subplot=1,label="RK41")     
     Plots.plot!(T,Y[:,2],subplot=2,label="RK41")

     println("RK42")
     T,Y=ode_rk42(fun_vdp,[t0 tf],y0,N[i])
     display("text/plain",[T Y])
     Plots.plot!(T,Y[:,1],subplot=1,label="RK42")     
     Plots.plot!(T,Y[:,2],subplot=2,label="RK42")
     # Gauss
     
     options[2] = 15
     # epsilon pour le test du point fixe
     # Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12 comme erreur absolue
     options[3] = 1e-6
     options[1] = N[i]
     display(options)
     T,Y,nphie,ifail=ode_gauss_pf(fun_vdp,[t0 tf],y0,options)
     display("text/plain",[T Y])
     display("text/plain",nphie)
     display("text/plain",ifail)
     Ifail[i] = length(findall(ifail==-1))
     display(plt)
end
display(Ifail)