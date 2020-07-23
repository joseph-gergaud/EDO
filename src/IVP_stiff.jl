#
# pb raides
#

using Plots 

include("functions/fun_stiff.jl")
include("solvers/ode_euler.jl")
include("solvers/ode_runge.jl")
include("solvers/ode_rk41.jl")
include("solvers/ode_rk42.jl")
include("solvers/ode_heun.jl")
include("solvers/ode_gauss_pf.jl")
include("solvers/ode_euler_pf.jl")

pause(text) = (print(stdout, text); read(stdin, 1); nothing)

y0 = [10]
t0 = 0
tf = 1.5

N0=[round(tf * 50 / 1.974 - 0.49) round(tf * 50 / 1.875 - 0.49) 50 100] 
N=N0
N = [30 40 80 100]

options = zeros(3)
for i=1:length(N)
    # Euler
    T,Y = ode_euler(fun_stiff,[t0 tf],y0,N[i])
    plt = Plots.plot(T,Y,color="blue",label="Euler",title="Résultats avec N="*string(N[i]))

    # Runge
    T,Y = ode_runge(fun_stiff,[t0 tf],y0,N[i])
    Plots.plot!(T,Y,color="cyan",label="Runge")

    # Heun
    T,Y = ode_heun(fun_stiff,[t0 tf],y0,N[i])
    Plots.plot!(T,Y,color="yellow",label="Heun")

    # RK4 classique
    T,Y = ode_rk41(fun_stiff,[t0 tf],y0,N[i])
    Plots.plot!(T,Y,color="green",label="RK41")

    # RK4 regle 3/8
    T,Y = ode_rk42(fun_stiff,[t0 tf],y0,N[i])
    Plots.plot!(T,Y,color="magenta",label="RK42")

    # Gauss
    options[2] = 40
    options[3] = 1e-6
    options[1] = N[i]
    T,Y,nphie,ifail = ode_gauss_pf(fun_stiff,[t0 tf],y0,options)
    Plots.plot!(T,Y,color="red", xlabel="t", ylabel="y(t)",label="Gauss")
    #ifail
    display(plt)
    pause("tapez entrée")
end

plt1 = Plots.plot(layout=(2,2),title="Gauss")
plt2 = Plots.plot(layout=(2,2),title="Euler implicite")
for i=1:length(N)

    # Gauss
    options[2] = 40
    options[3] = 1e-6
    options[1] = N[i]

    T,Y,nphie,ifail = ode_gauss_pf(fun_stiff,[t0 tf],y0,options)
    display(ifail)
    Plots.plot!(plt1,T,Y,color="red",subplot=i,label="N="*string(N[i]))

    T,Y,nphie,ifail = ode_euler_pf(fun_stiff,[t0 tf],y0,options)
    display(ifail)
    Plots.plot!(plt2,T,Y,color="blue",subplot=i,label="N="*string(N[i]))
    # print("-depsc",fich)
end
display(plt1)
pause("tapez entrée")
display(plt2)
#=
for i=1:4
    figure(i)
    fich=["figstiff" int2str(i)]
    print("-depsc",fich)
    unix(["epstopdf figstiff" int2str(i) ".eps"])
end
=#
