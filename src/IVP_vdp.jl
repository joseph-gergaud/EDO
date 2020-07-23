#
# Intégration de l'equation differentiel de l'equation de Van der Pol
# ref: Hairer
#
using Plots
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
# fonction pause
pause(text) = (print(stdout, text); read(stdin, 1); nothing)

y0 = [2.008619860874843136; 0]
t0 = 0
tf = 6.663286859323130189
N0 = [120:60:1080; 1200:600:10800 ]# 12000:12000:72000] # multiple de 12 pour avoir des nombres entier si on divise
options = zeros(3)
# par 3 ou 4
#
# Solutions y_1 et y_2 et plan de phase pour RKE
# ----------------------------------------------
N = 25
pyplot()
plt = Plots.plot(layout=(3))

println("Euler")
println("-----")
T,Y = ode_euler(fun_vdp,[t0 tf],y0,N)
println("[T Y]=")
display("text/plain",[T Y])
plot_sol(plt,T,Y,"blue","euler")

println("Runge")
println("-----")
T,Y = ode_runge(fun_vdp,[t0 tf],y0,N)
println("[T Y]=")
display("text/plain",[T Y])
plot_sol(plt,T,Y,"olive","runge")

println("Heun")
println("----")
T,Y = ode_heun(fun_vdp,[t0 tf],y0,N)
println("[T Y]=")
display("text/plain",[T Y])
plot_sol(plt,T,Y,"magenta","heun")

display("text/plain","rk41")
println("----")
T,Y = ode_rk41(fun_vdp,[t0 tf],y0,N)
println("[T Y]=")
display("text/plain",[T Y])
plot_sol(plt,T,Y,"green","rk4_1")

println("rk42")
println("-----")
T,Y = ode_rk42(fun_vdp,[t0 tf],y0,N)
println("[T Y]=")
display("text/plain",[T Y])
plot_sol(plt,T,Y,"cyan","rk4_2")

println("Gauss, point fixe")
println("-----------------")
options[1] = N
# options[2] = nbre max iteration point fixe
options[2] = 15
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12 comme erreur absolue
options[3] = 1e-6
println("[N, nb_itmax, f_eps]")
display("text/plain",options)
T,Y,nfe,ifail = ode_gauss_pf(fun_vdp,[t0 tf],y0,options)
println("[T Y]=")
display("text/plain",[T Y])
display(nfe)
display(ifail')
plot_sol(plt,T,Y,"red","gauss_pf")

println("Gauss, Newton")
println("-------------")
println("[N, nb_itmax, f_eps]")
display("text/plain",options')
T,Y,nfe,ndphie,ifail = ode_gauss_newton(fun_vdp,d_fun_vdp,[t0 tf],y0,options)
println("[T Y]=")
display("text/plain",[T Y])
println("nfe=")
display("text/plain",nfe)
println("ndphie=")
display("text/plain",ndphie)
println("ifail=")
display("text/plain",ifail')
plot_sol(plt,T,Y,"red","Gauss_Newton")

#print -depsc fig_solutions_vdp
#format short
#diary

pause("tapez entrée pour voir le graphique des ordres")

# Courbes d"ordre
# ---------------

N = N0
# Euler
closeall()
err1 = zeros(length(N))
err2 = zeros(length(N))
Nfe = zeros(length(N))
Ifail = zeros(length(N))

plt=Plots.plot(layout=(1,2))
for i=1:length(N)
    T,Y = ode_euler(fun_vdp,[t0 tf],y0,N[i])
    n = length(T)
    err1[i] = abs(y0[1]-Y[n,1])
    err2[i] = abs(y0[2]-Y[n,2])
end
#
v = [log10(N[1]),log10(err1[1]),log10(N[end]),log10(err1[end])]
Plots.plot!(log10.(N),log10.(err1), xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="Euler")
Plots.plot!(log10.(N),log10.(err2), xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="Euler")
i = rand(1:length(N))
annotate!([(log10.(N[i]),log10.(err1[i]), Plots.text("Euler", 10, :center))],subplot=1)
annotate!([(log10.(N[i]),log10.(err2[i]), Plots.text("Euler", 10,:center))],subplot=2)

#text(log10.(N[i]),log10.(err2[i] ) ,"Euler") 
#
# Runge
s = 2
N = N0/s

for i=1:length(N)
    T,Y = ode_runge(fun_vdp,[t0 tf],y0,N[i])
    n = length(T)
    err1[i] = abs(y0[1]-Y[n,1])
    err2[i] = abs(y0[2]-Y[n,2])
end

Plots.plot!(log10.(s*N),log10.(err1),color="olive", xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="Runge")
Plots.plot!(log10.(s*N),log10.(err2),color="olive", xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="Runge")
i = rand(1:length(N))
annotate!([(log10.(s*N[i]),log10.(err1[i]), Plots.text("Runge", 10,:olive, :center))],subplot=1)
annotate!([(log10.(s*N[i]),log10.(err2[i]), Plots.text("Runge", 10,:olive, :center))],subplot=2)
# Heun
s = 3
N = N0/s

for i=1:length(N)
    T,Y = ode_heun(fun_vdp,[t0 tf],y0,N[i])
    n = length(T)
    err1[i]=abs(y0[1]-Y[n,1])
    err2[i]=abs(y0[2]-Y[n,2])
end
   
Plots.plot!(log10.(s*N),log10.(err1),color="magenta", xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="Heun")
Plots.plot!(log10.(s*N),log10.(err2),color="magenta", xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="Heun")
i = rand(1:length(N))
annotate!([(log10.(s*N[i]),log10.(err1[i]), Plots.text("Heun", 10,:magenta, :center))],subplot=1)
annotate!([(log10.(s*N[i]),log10.(err2[i]), Plots.text("Heun", 10,:magenta, :center))],subplot=2)

# RK4 classique
s = 4
N = N0/s

for i=1:length(N)
    T,Y = ode_rk41(fun_vdp,[t0 tf],y0,N[i])
    n = length(T)
    err1[i] = abs(y0[1]-Y[n,1])
    err2[i] = abs(y0[2]-Y[n,2])
end

Plots.plot!(log10.(s*N),log10.(err1),color="green", xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="RK41")
Plots.plot!(log10.(s*N),log10.(err2),color="green", xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="RK41")
i = rand(1:length(N))
annotate!([(log10.(s*N[i]),log10.(err1[i]), Plots.text("RK41", 10,:green, :center))],subplot=1)
annotate!([(log10.(s*N[i]),log10.(err2[i]), Plots.text("RK41", 10,:green, :center))],subplot=2)

# RK4 regle 3/8
s = 4
N = N0/s

for i=1:length(N)
    T,Y = ode_rk42(fun_vdp,[t0 tf],y0,N[i])
    n = length(T)
    err1[i]=abs(y0[1]-Y[n,1])
    err2[i]=abs(y0[2]-Y[n,2])
end

Plots.plot!(log10.(s*N),log10.(err1),color="cyan", xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="RK42")
Plots.plot!(log10.(s*N), log10.(err2),color="cyan", xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="RK42")
i = rand(1:length(N))
annotate!([(log10.(s*N[i]),log10.(err1[i]), Plots.text("RK42", 10,:cyan, :center))],subplot=1)
annotate!([(log10.(s*N[i]),log10.(err2[i]), Plots.text("RK42", 10,:cyan, :center))],subplot=2)

# print -depsc fig1_ordre
# print -depsc fig2_ordre
#
# Gauss + point fixe
# -----
s=4
N=N0/s
# options[2] = nbre max iteration point fixe
options[2] = 15
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12 comme erreur absolue
options[3] = 1e-12
Nfe = zeros(length(N))
Ifail = zeros(length(N))
for i=1:length(N)
    options[1] = N[i]
    T,Y,nfe,ifail = ode_gauss_pf(fun_vdp,[t0 tf],y0,options)    
    Ifail[i] = length(findall(ifail==-1)) 
    n = length(T)
    err1[i] = abs(y0[1]-Y[n,1])
    err2[i] = abs(y0[2]-Y[n,2])
    Nfe[i] = nfe
end

#print -depsc fig_ordre_Gauss1
# vraie courbe
# ------------
Plots.plot!(log10.(Nfe),log10.(err1),color="red",xlabel="log_{10}(fe)",ylabel="log_{10}(erreur pour y_1)",subplot=1,label="Gauss pf")
Plots.plot!(log10.(Nfe), log10.(err2),color= "red", xlabel="log_{10}(nphi)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="Gauss pf")
i = rand(1:length(N))
annotate!([(log10.(Nfe[i]),log10.(err1[i]), Plots.text("Gauss pf", 10,:red, :center))],subplot=1)
annotate!([(log10.(Nfe[i]),log10.(err2[i]), Plots.text("Gauss pf", 10,:red, :center))],subplot=2)
##legend("Euler", "Runge","Heun","RK4_1","RK4_2","gauss pf feps=1.e-12","Location","SouthWest")
#print -depsc fig_ordre_Gauss2
display(plt)
pause("tapez entrée pour voir les vraies Courbes")
plt = Plots.plot(layout=(1,2))

# vraie courbe
# ------------
Plots.plot!(log10.(Nfe),log10.(err1),color="red", xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="Gauss pf feps=1e-12")
Plots.plot!(log10.(Nfe),log10.(err2),color="red", xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="Gauss pf feps=1e-12")
i = rand(1:length(N))
annotate!([(log10.(Nfe[i]),log10.(err1[i]), Plots.text("Gauss pf feps=1e-12", 10,:red, :center))],subplot=1)
annotate!([(log10.(Nfe[i]),log10.(err2[i]), Plots.text("Gauss pf feps=1e-12", 10,:red, :center))],subplot=2)
#
# Gauss feps=1.e-6
s = 4
N = N0/s
# options[2] = nbre max iteration point fixe
options[2] = 15
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12 comme erreur absolue
options[3] = 1e-6

for i=1:length(N)
    options[1] = N[i]
    T,Y,nfe,ifail = ode_gauss_pf(fun_vdp,[t0 tf],y0,options)    
    Ifail[i] = length(findall(ifail==-1)) 
    n = length(T)
    err1[i] = abs(y0[1]-Y[n,1])
    err2[i] = abs(y0[2]-Y[n,2])
    Nfe[i] = nfe
end
# 
# vraie courbe
# ------------
Plots.plot!(log10.(Nfe),log10.(err1), xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="Gauss pf feps=1e-6")
Plots.plot!(log10.(Nfe),log10.(err2), xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="Gauss pf feps=1e-6")
i = rand(1:length(N))
annotate!([(log10.(Nfe[i]),log10.(err1[i]), Plots.text("Gauss pf feps=1e-6", 10,:center))],subplot=1)
annotate!([(log10.(Nfe[i]),log10.(err2[i]), Plots.text("Gauss pf feps=1e-6", 10,:center))],subplot=2)
#
# Gauss fpitermax=2
s = 4
N = N0/s
# options[2] = nbre max iteration point fixe
options[2] = 2
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12 comme erreur absolue
options[3] = 1e-12

for i=1:length(N)
    options[1] = N[i]
    T,Y,nfe,ifail = ode_gauss_pf(fun_vdp,[t0 tf],y0,options)    
    Ifail[i] = length(findall(ifail==-1)) 
    n = length(T)
    err1[i] = abs(y0[1]-Y[n,1])
    err2[i] = abs(y0[2]-Y[n,2])
    Nfe[i] = nfe
end
# vraie courbe
# ------------
Plots.plot!(log10.(Nfe),log10.(err1),color="green",xlabel="log_{10}(nfe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="Gauss pf fpitermax=2")
Plots.plot!(log10.(Nfe),log10.(err2),color="green",xlabel="log_{10}(fnphi)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="Gauss pf fpitermax=2")
i = rand(1:length(N))
annotate!([(log10.(Nfe[i]),log10.(err1[i]), Plots.text("Gauss pf fpitermax=2", 10,:green, :center))],subplot=1)
annotate!([(log10.(Nfe[i]),log10.(err2[i]), Plots.text("Gauss pf fpitermax=2", 10,:green, :center))],subplot=2)
#legend("gauss pf feps=1.e-12","gauss pf feps=1.e-6","gauss pf fpitermax=2","Location","SouthWest")
#print -depsc fig_ordre_Gauss3
#
# Gauss + Newton
# ---------------
#=
s = 2
N = N0/(2*s)
# options[2] = nbre max iteration point fixe
options[2] = 15
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12
# comme erreur absolue
options[3] = 1e-12

for i=1:length(N)
    options[1] = N[i]
    T,Y,nfe,ndphie,ifail = ode_gauss_newton(fun_vdp,d_fun_vdp,[t0 tf],y0,options)    
    Ifail[i] = length(findall(ifail==-1)) 
    n = length(T)
    err1[i] = abs(y0[1]-Y[n,1])
    err2[i] = abs(y0[2]-Y[n,2])
    Nfe[i] = nfe + n * ndphie
end

# vraie courbe
# ------------
Plots.plot!(log10.(Nfe),log10.(err1),color="blue",ylabel="log_{10}(erreur pour y_1)",subplot=1,label="Gauss Newton feps=1e-12")
Plots.plot!(log10.(Nfe),log10.(err2),color="blue",ylabel="log_{10}(erreur pour y_2)",subplot=2,label="Gauss Newton feps=1e-12")
i = rand(1:length(N))
annotate!([(log10.(Nfe[i]),log10.(err1[i]), Plots.text("Gauss Newton feps=1e-12", 10,:blue, :center))],subplot=1)
annotate!([(log10.(Nfe[i]),log10.(err2[i]), Plots.text("Gauss Newton feps=1e-12", 10,:blue, :center))],subplot=2)
#legend("Euler", "Runge","Heun","RK4_1","RK4_2","gauss pf feps=1.e-12","Location","SouthWest")
=#
display(plt)