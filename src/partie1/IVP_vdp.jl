using Plots
using LinearAlgebra
"""
~gergaud/ENS/Control/ODE/Matlab/IVP_vdp.m

Auteurs:  Joseph GERGAUD
Date:     nov. 2005
Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
        2, rue Camichel 31071 Toulouse FRANCE
Email:    gergaudenseeiht.fr
***************************************************************************

Int�gration de l'equation differentiel de l'equation de Van der Pol
ref: Hairer


Initialisations
---------------
""" 

include("ode_gauss_v3.jl")
include("ode_gauss_v2.jl")
include("ode_gauss.jl")
include("edo_euler.jl")
include("edo_runge.jl")
include("edo_rk41.jl")
include("edo_rk42.jl")
include("edo_heun.jl")
include("plot_sol.jl")
include("d_phi_vdp.jl")
include("phi_vdp.jl")

#nom_fich = "IVP_vdp_N_10.txt"
#if exist(nom_fich,"file")
#    delete(nom_fich)
#end
#format long
#diary(nom_fich)
#
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
N = 10
#println("Euler")
#println("-----")
pyplot()
plt=Plots.plot(layout=(1,3))
T,Y = ode_euler(phi_vdp,[t0 tf],y0,N)
#println("[T Y]=")
#println([T Y])
plot_sol(plt,T,Y,"blue","euler")
#println("Runge")
#println("-----")
T,Y = ode_runge(phi_vdp,[t0 tf],y0,N)
#println("[T Y]=")
#println([T Y])
plot_sol(plt,T,Y,"olive","runge")
#println("Heun")
#println("----")
T,Y = ode_heun(phi_vdp,[t0 tf],y0,N)
#println("[T Y]=")
#println([T Y])
plot_sol(plt,T,Y,"magenta","heun")
#println("rk41")
#println("----")
T,Y = ode_rk41(phi_vdp,[t0 tf],y0,N)
#println("[T Y]=")
#println([T Y])
plot_sol(plt,T,Y,"green","rk4_1")
#println("rk42")
#println("-----")
T,Y = ode_rk42(phi_vdp,[t0 tf],y0,N)
#println("[T Y]=")
#println([T Y])
plot_sol(plt,T,Y,"cyan","rk4_2")
#println("Gauss, point fixe")
#println("-----------------")
options[1] = N
# options[2] = nbre max iteration point fixe
options[2] = 15
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12 comme erreur absolue
options[3] = 1e-6
#println("[N, nb_itmax, f_eps]")
#println(options)
T,Y,nphie,ifail = ode_gauss_v2(phi_vdp,[t0 tf],y0,options)
#println("[T Y]=")
#println([T Y])
#nphie
#ifail
plot_sol(plt,T,Y,"red","gauss_v2")

#println("Gauss, Newton")
#println("-------------")
#println("[N, nb_itmax, f_eps]")
#println(options)
T,Y,nphie,ndphie,ifail = ode_gauss_v3(phi_vdp,d_phi_vdp,[t0 tf],y0,options)
#println("[T Y]=")
#println([T Y])
#nphie
#ndphie
#ifail
plot_sol(plt,T,Y,"red","Gauss_Newton")

#print -depsc fig_solutions_vdp
#format short
#diary
#
pause("tapez entrée pour voir le graphique des ordres")

# Courbes d"ordre
# ---------------


N = N0
# Euler
err1 = zeros(length(N))
err2 = zeros(length(N))
plt=Plots.plot(layout=(1,2))
for i=1:length(N)
    T,Y = ode_euler(phi_vdp,[t0 tf],y0,N[i])
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
#clear T Y err1 err2
s = 2
N = N0/s
err1 = zeros(length(N))
err2 = zeros(length(N))

for i=1:length(N)
    T,Y = ode_runge(phi_vdp,[t0 tf],y0,N[i])
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
#clear T Y err1 err2
s = 3
N = N0/s
err1 = zeros(length(N))
err2 = zeros(length(N))

for i=1:length(N)
    T,Y = ode_heun(phi_vdp,[t0 tf],y0,N[i])
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
##clear T Y err1 err2
s = 4
N = N0/s
err1 = zeros(length(N))
err2 = zeros(length(N))

for i=1:length(N)
    T,Y = ode_rk41(phi_vdp,[t0 tf],y0,N[i])
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
err1 = zeros(length(N))
err2 = zeros(length(N))

for i=1:length(N)
    T,Y = ode_rk42(phi_vdp,[t0 tf],y0,N[i])
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
Nphie = []
s=4
N=N0/s
# options[2] = nbre max iteration point fixe
options[2] = 15
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12 comme erreur absolue
options[3] = 1e-12
Ifail = []
err1 = zeros(length(N))
err2 = zeros(length(N))
Nphie = zeros(length(N))
Ifail = zeros(length(N))
for i=1:length(N)
    options[1] = N[i]
    T,Y,nphie,ifail = ode_gauss_v2(phi_vdp,[t0 tf],y0,options)    
    Ifail[i] = length(findall(ifail==-1)) 
    n = length(T)
    err1[i] = abs(y0[1]-Y[n,1])
    err2[i] = abs(y0[2]-Y[n,2])
    Nphie[i] = nphie
end

#print -depsc fig_ordre_Gauss1
# vraie courbe
# ------------
Plots.plot!(log10.(Nphie),log10.(err1),color="red",xlabel="log_{10}(fe)",ylabel="log_{10}(erreur pour y_1)",subplot=1,label="Gauss")
Plots.plot!(log10.(Nphie), log10.(err2),color= "red", xlabel="log_{10}(nphi)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="Gauss")
i = rand(1:length(N))
annotate!([(log10.(Nphie[i]),log10.(err1[i]), Plots.text("Gauss", 10,:red, :center))],subplot=1)
annotate!([(log10.(Nphie[i]),log10.(err2[i]), Plots.text("Gauss", 10,:red, :center))],subplot=2)
##legend("Euler", "Runge","Heun","RK4_1","RK4_2","gauss v2 feps=1.e-12","Location","SouthWest")
#print -depsc fig_ordre_Gauss2
display(plt)
pause("tapez entrée pour voir les vraies Courbes")
plt = Plots.plot(layout=(1,2))
# vraie courbe
# ------------
Plots.plot!(log10.(Nphie),log10.(err1),color="red", xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="Gauss v2 feps=1e-12")
Plots.plot!(log10.(Nphie),log10.(err2),color="red", xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="Gauss v2 feps=1e-12")
i = rand(1:length(N))
annotate!([(log10.(Nphie[i]),log10.(err1[i]), Plots.text("Gauss v2 feps=1e-12", 10,:red, :center))],subplot=1)
annotate!([(log10.(Nphie[i]),log10.(err2[i]), Plots.text("Gauss v2 feps=1e-12", 10,:red, :center))],subplot=2)
#
# Gauss feps=1.e-6
Nphie = []
s = 4
N = N0/s
# options[2] = nbre max iteration point fixe
options[2] = 15
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12 comme erreur absolue
options[3] = 1e-6
Ifail = []
err1 = zeros(length(N))
err2 = zeros(length(N))
Nphie = zeros(length(N))
Ifail = zeros(length(N))
for i=1:length(N)
    options[1] = N[i]
    T,Y,nphie,ifail = ode_gauss_v2(phi_vdp,[t0 tf],y0,options)    
    Ifail[i] = length(findall(ifail==-1)) 
    n = length(T)
    err1[i] = abs(y0[1]-Y[n,1])
    err2[i] = abs(y0[2]-Y[n,2])
    Nphie[i] = nphie
end
# 
# vraie courbe
# ------------
Plots.plot!(log10.(Nphie),log10.(err1), xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="Gauss v2 feps=1e-6")
Plots.plot!(log10.(Nphie),log10.(err2), xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="Gauss v2 feps=1e-6")
i = rand(1:length(N))
annotate!([(log10.(Nphie[i]),log10.(err1[i]), Plots.text("Gauss v2 feps=1e-6", 10,:center))],subplot=1)
annotate!([(log10.(Nphie[i]),log10.(err2[i]), Plots.text("Gauss v2 feps=1e-6", 10,:center))],subplot=2)
#
# Gauss fpitermax=2
#clear T Y err1 err2
Nphie = []
s = 4
N = N0/s
# options[2] = nbre max iteration point fixe
options[2] = 2
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12 comme erreur absolue
options[3] = 1e-12
Ifail = []
err1 = zeros(length(N))
err2 = zeros(length(N))
Nphie = zeros(length(N))
Ifail = zeros(length(N))
for i=1:length(N)
    options[1] = N[i]
    T,Y,nphie,ifail = ode_gauss_v2(phi_vdp,[t0 tf],y0,options)    
    Ifail[i] = length(findall(ifail==-1)) 
    n = length(T)
    err1[i] = abs(y0[1]-Y[n,1])
    err2[i] = abs(y0[2]-Y[n,2])
    Nphie[i] = nphie
end
# vraie courbe
# ------------
Plots.plot!(log10.(Nphie),log10.(err1),color="green",xlabel="log_{10}(nphie)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="Gauss v2 fpitermax=2")
Plots.plot!(log10.(Nphie),log10.(err2),color="green",xlabel="log_{10}(fnphi)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="Gauss v2 fpitermax=2")
i = rand(1:length(N))
annotate!([(log10.(Nphie[i]),log10.(err1[i]), Plots.text("Gauss v2 fpitermax=2", 10,:green, :center))],subplot=1)
annotate!([(log10.(Nphie[i]),log10.(err2[i]), Plots.text("Gauss v2 fpitermax=2", 10,:green, :center))],subplot=2)
#legend("gauss v2 feps=1.e-12","gauss v2 feps=1.e-6","gauss v2 fpitermax=2","Location","SouthWest")
#print -depsc fig_ordre_Gauss3
#
# Gauss + Newton
# ---------------
Nphie = []
s=2
N=N0/(2*s)
# options[2] = nbre max iteration point fixe
options[2] = 15
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12
# comme erreur absolue
options[3] = 1e-12
Ifail = []
err1 = zeros(length(N))
err2 = zeros(length(N))
Nphie = zeros(length(N))
Ifail = zeros(length(N))
for i=1:length(N)
    options[1] = N[i]
    T,Y,nphie,ndphie,ifail = ode_gauss_v3(phi_vdp,d_phi_vdp,[t0 tf],y0,options)    
    Ifail[i] = length(findall(ifail==-1)) 
    n = length(T)
    err1[i] = abs(y0[1]-Y[n,1])
    err2[i] = abs(y0[2]-Y[n,2])
    Nphie[i] = nphie + n * ndphie
end

# vraie courbe
# ------------
#Plots.plot!(log10.(Nphie),log10.(err1),color="blue",ylabel="log_{10}(erreur pour y_1)",subplot=1,label="Gauss v3 feps=1e-12")
#Plots.plot!(log10.(Nphie),log10.(err2),color="blue",ylabel="log_{10}(erreur pour y_2)",subplot=2,label="Gauss v3 feps=1e-12")
i = rand(1:length(N))
#annotate!([(log10.(Nphie[i]),log10.(err1[i]), Plots.text("Gauss v3 feps=1e-12", 10,:blue, :center))],subplot=1)
#annotate!([(log10.(Nphie[i]),log10.(err2[i]), Plots.text("Gauss v3 feps=1e-12", 10,:blue, :center))],subplot=2)
##legend("Euler", "Runge","Heun","RK4_1","RK4_2","gauss v2 feps=1.e-12","Location","SouthWest")
display(plt)