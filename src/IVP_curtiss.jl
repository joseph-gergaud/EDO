#
# Intégration de l'equation differentiel de l'équation Curtiss et Hirschfelder
# ref: Hairer page 3 tome 2; pb raide
#

using Plots
using LinearAlgebra

include("functions/fun_curtiss.jl")
include("functions/fun_toto.jl")
include("solvers/ode_euler_pf.jl")
include("solvers/ode_euler.jl")
include("solvers/ode_runge.jl")
include("solvers/ode_heun.jl")
include("solvers/ode_rk41.jl")
include("solvers/ode_rk42.jl")
include("solvers/ode_gauss_pf.jl")

pause(text) = (print(stdout, text); read(stdin, 1); nothing)

pyplot()
closeall()

y0 = [0]
t0 = 0
tf = 1.5
N0 = [floor(tf*50/1.974), floor(tf*50/1.875), 50, 100]
N = N0
options = zeros(3,1)


for i=1:length(N)
    #figure(i)

    pause("tapez Entrée pour voir les schémas (N = "*string(N[i])*" )")
    plt = Plots.plot(title="N = "*string(N[i]))
    # Euler
    T,Y = ode_euler(fun_curtiss,[t0 tf],y0,N[i])
    Plots.plot!(T,Y,LineWidth=2,color=:blue, label = "Euler")
    
    # Runge
    T,Y = ode_runge(fun_curtiss,[t0 tf],y0,N[i])
    Plots.plot!(T,Y,LineWidth = 2,color =:black,label = "Runge")
   
    # Heun
    T,Y = ode_heun(fun_curtiss,[t0 tf],y0,N[i])
    Plots.plot!(T,Y,LineWidth = 2,color=:yellow,label = "Heun")

    # RK4 classique
    T,Y = ode_rk41(fun_curtiss,[t0 tf],y0,N[i])
    Plots.plot!(T,Y,LineWidth=2,label = "RK41")

    # RK4 regle 3/8
    T,Y = ode_rk42(fun_curtiss,[t0 tf],y0,N[i])
    Plots.plot!(T,Y,LineWidth = 2,label = "RK42" )

    # Gauss
    options[2] = 40
    options[3] = 1e-6
    options[1] = N[i]
    T,Y,nphie,ifail = ode_gauss_pf(fun_curtiss,[t0 tf],y0,options)
    println(ifail)
    Plots.plot!(T,Y,LineWidth = 2,color=:red,xlabel = "t",ylabel = "y(t)" , label = "Gauss")
    display(plt)
end

# Gauss
plt1 = Plots.plot(layout = (2,2))
# Euler implicite
plt2 = Plots.plot(layout = (2,2))

for i=1:length(N)
	options[2] = 40
	options[3] = 1e-6
	options[1] = N[i]
	# Gauss	
	T,Y,nphie,ifail = ode_gauss_pf(fun_curtiss,[t0 tf],y0,options)
	println(ifail)
	Plots.plot!(plt1[i],T,Y,color=:red,LineWidth=2,label="Gauss,N = "*string(N[i]))
	# Euler implicite
	T,Y,nphie,ifail = ode_euler_pf(fun_curtiss,[t0 tf],y0,options)
	println(ifail)
	Plots.plot!(plt2[i],T,Y,color=:green,LineWidth=2,label="Euler implicite,N = "*string(N[i]))
end

pause("tapez Entrée pour voir les vraies Figures avec le schéma de Gauss")
display(plt1)
pause("tapez Entrée pour voir les vraies Figures avec le schéma d'euler implicite")
display(plt2)
#=
for i=1:6
    figure(i)
    fich=["figcurtiss" int2str(i)]
    print("-depsc",fich)
    unix(["epstopdf figcurtiss' int2str(i) '.eps"])
end
=#
