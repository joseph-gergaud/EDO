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
include("plot_sol.jl")
include("d_phi_vdp.jl")
include("phi_vdp.jl")
include("edo_euler.jl")
include("edo_runge.jl")
include("edo_rk4.jl")
include("edo_heun.jl")
#clear all close all
#nom_fich = "IVP_vdp_N_10.txt"
#if exist(nom_fich,"file")
#    delete(nom_fich)
#end
#format long
#diary(nom_fich)
#


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
#figure(4)
#println("Euler")
#println("-----")
T,Y = ode_euler(phi_vdp,[t0 tf],y0,N)
#println("[T Y]=")
#println([T Y])
plot_sol(T,Y,"blue",["euler","euler","euler"])
#println("Runge")
#println("-----")
T,Y = ode_runge(phi_vdp,[t0 tf],y0,N)
#println("[T Y]=")
#println([T Y])
plot_sol(T,Y,"yellow",["runge","runge","runge"])
#println("Heun")
#println("----")
T,Y = ode_heun(phi_vdp,[t0 tf],y0,N)
#println("[T Y]=")
#println([T Y])
plot_sol(T,Y,"magenta",["heun","heun","heun"])
#println("rk41")
#println("----")
T,Y = ode_rk4(phi_vdp,[t0 tf],y0,N)
#println("[T Y]=")
#println([T Y])
plot_sol(T,Y,"green",["rk41","rk41","rk41"])
#println("rk42")
#println("-----")
T,Y = ode_rk4(phi_vdp,[t0 tf],y0,N)
#println("[T Y]=")
#println([T Y])
plot_sol(T,Y,"cyan",["rk42","rk42","rk42"])
#println("Gauss, point fixe")
#println("-----------------")
options[1] = N
# options(2) = nbre max iteration point fixe
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
plot_sol(T,Y,"red",["gauss_v2","gauss_v2","gauss_v2"])

#println("Gauss, Newton")
#println("-------------")
#println("[N, nb_itmax, f_eps]")
#println(options)
T,Y,nphie,ndphie,ifail = ode_gauss_v3(phi_vdp,d_phi_vdp,[t0 tf],y0,options)
#println("[T Y]=")
#println([T Y])
nphie
ndphie
ifail
plot_sol(T,Y,"red",["gauss_v3","gauss_v3","gauss_v3"])

#legend("Euler", "Runge","Heun","RK4_1","RK4_2","gauss v2","gauss_Newton","Location","SouthEast")


#print -depsc fig_solutions_vdp
#format short
#diary
#
#pause

# Courbes d"ordre
# ---------------


N=N0
# Euler
for i=1:length(N),
    T,Y = ode_euler(phi_vdp,[t0 tf],y0,N(i))
    n=length(T)
    err1(i)=abs(y0(1)-Y(n,1))
    err2(i)=abs(y0(2)-Y(n,2))
end
figure(1)
subplot(1,2,1)
loglog(N,err1)
i=round(length(N)/2)
text(N(i),err1(i),"Euler")
xlabel("fe")
ylabel("erreur pour y_1")     
subplot(1,2,2)
loglog(N,err2)
text(N(i),err2(i),"Euler") 
xlabel("fe")
ylabel("erreur pour y_2")    
#
figure(2)
subplot(1,2,1)
Plots.plot(log10(N),log10(err1))
i=round(length(N)/2)
text(log10(N(i)),log10(err1(i)),"Euler")
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_1)")    
subplot(1,2,2)
Plots.plot(log10(N),log10(err2))
text(log10(N(i)),log10(err2(i)),"Euler") 
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_2)")  
#
# Runge
#clear T Y err1 err2
s=2
N=N0/s
for i=1:length(N),
T,Y = ode_runge(phi_vdp,[t0 tf],y0,N(i))
n=length(T)
err1(i)=abs(y0(1)-Y(n,1))
err2(i)=abs(y0(2)-Y(n,2))
end
figure(1)
subplot(1,2,1)
#hold on
loglog(s*N,err1)
i=round(length(N)/2)
#text(s*N(i),err1(i),"Runge")    
subplot(1,2,2)
#hold on
loglog(s*N,err2)
text(s*N(i),err2(i),"Runge")  
figure(2)
subplot(1,2,1)
#hold on
    Plots.plot(log10(s*N),log10(err1),"k")
i=round(length(N)/2)
#text(log10(s*N(i)),log10(err1(i)),"Runge")
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_1)")     
subplot(1,2,2)
#hold on
    Plots.plot(log10(s*N),log10(err2),"k")
#     text(log10(s*N(i)),log10(err2(i)),"Runge") 
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_2)")    
# Heun
#clear T Y err1 err2
s=3
N=N0/s
for i=1:length(N),
    T,Y = ode_heun(phi_vdp,[t0 tf],y0,N(i))
    n=length(T)
    err1(i)=abs(y0(1)-Y(n,1))
    err2(i)=abs(y0(2)-Y(n,2))
end
figure(1)
subplot(1,2,1)
#hold on
loglog(s*N,err1)
i=round(length(N)/2)
text(s*N(i),err1(i),"Heun")    
subplot(1,2,2)
#hold on
loglog(s*N,err2)
text(s*N(i),err2(i),"Heun")    
figure(2)
subplot(1,2,1)
#hold on
Plots.plot(log10(s*N),log10(err1),"m")
i=round(length(N)/2)
#text(log10(s*N(i)),log10(err1(i)),"Heun")
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_1)")     
subplot(1,2,2)
#hold on
Plots.plot(log10(s*N),log10(err2),"m")
#     text(log10(s*N(i)),log10(err2(i)),"Heun") 
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_2)")    

# RK4 classique
##clear T Y err1 err2
s=4
N=N0/s
for i=1:length(N),
    T,Y = ode_rk4(phi_vdp,[t0 tf],y0,N(i))
    n=length(T)
    err1(i)=abs(y0(1)-Y(n,1))
    err2(i)=abs(y0(2)-Y(n,2))
end
figure(1)
subplot(1,2,1)
#hold on
loglog(s*N,err1)
i=round(length(N)/2)
text(s*N(i),err1(i),"RK4")    
subplot(1,2,2)
#hold on
loglog(s*N,err2)
text(s*N(i),err2(i),"RK4") 
figure(2)
subplot(1,2,1)
#hold on
Plots.plot(log10(s*N),log10(err1),"g")
i=round(length(N)/2)
#text(log10(s*N(i)),log10(err1(i)),"Rk4")
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_1)")     
subplot(1,2,2)
#hold on
Plots.plot(log10(s*N),log10(err2),"g")
#     text(log10(s*N(i)),log10(err2(i)),"Rk4") 
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_2)")    

# RK4 regle 3/8
#clear T Y err1 err2
s=4
N=N0/s
for i=1:length(N),
    T,Y = ode_rk4(phi_vdp,[t0 tf],y0,N(i))
    n=length(T)
    err1(i)=abs(y0(1)-Y(n,1))
    err2(i)=abs(y0(2)-Y(n,2))
end
figure(1)
subplot(1,2,1)
#hold on
loglog(s*N,err1)    
subplot(1,2,2)
#hold on
loglog(s*N,err2)
figure(2)
subplot(1,2,1)
#hold on
Plots.plot(log10(s*N),log10(err1),"c")
i=round(length(N)/2)
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_1)")     
subplot(1,2,2)
#hold on
Plots.plot(log10(s*N),log10(err2),"c")
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_2)")    
figure(1)
# print -depsc fig1_ordre
figure(2)
# print -depsc fig2_ordre
#
# Gauss + point fixe
# -----
##clear T Y err1 err2
Nphie = []
s=4
N=N0/s
# options(2) = nbre max iteration point fixe
options(2) = 15
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12 comme erreur absolue
options(3) = 1e-12
Ifail = []
for i=1:length(N),
    options(1) = N(i)
    [T,Y,nphie,ifail] = ode_gauss_v2(phi_vdp,[t0 tf],y0,options)
    Ifail = [Ifail length(find(ifail==-1))] 
    n=length(T)
    err1(i)=abs(y0(1)-Y(n,1))
    err2(i)=abs(y0(2)-Y(n,2))
    Nphie(i)=nphie
end
# Inutile de trie car deja en ordre
figure(1) # pour verifier l'ordre
 # ---------------------
subplot(1,2,1)
#hold on
loglog(s*N,err1,"red")    
subplot(1,2,2)
#hold on
loglog(s*N,err2,"red")
print -depsc fig_ordre_Gauss1
figure(2) 
# vraie courbe
# ------------
subplot(1,2,1)
Plots.plot(log10(Nphie),log10(err1),"red")
i=round(length(N)/2)
#text(log10(Nphie(i)),log10(err1(i)),"\color{red}Gauss")
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_1)")     
v=axis
axis([v(1) v(2) v(3)-4 v(4)])
subplot(1,2,2)
Plots.plot(log10(Nphie),log10(err2),"red")
i=round(length(N)/2)
#text(log10(Nphie(i)),log10(err2(i)),"\color{red}Gauss")
xlabel("log_{10}(nphi)")
ylabel("log_{10}(erreur pour y_2)") 
axis([v(1) v(2) v(3)-4 v(4)])
legend("Euler", "Runge","Heun","RK4_1","RK4_2","gauss v2 feps=1.e-12","Location","SouthWest")
print -depsc fig_ordre_Gauss2
figure(3) 
# vraie courbe
# ------------
subplot(1,2,1)
Plots.plot(log10(Nphie),log10(err1),"red")
i=round(length(N)/2)
#text(log10(Nphie(i)),log10(err1(i)),"\color{red}Gauss")
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_1)")     
subplot(1,2,2)
Plots.plot(log10(Nphie),log10(err2),"red")
i=round(length(N)/2)
#text(log10(Nphie(i)),log10(err2(i)),"\color{red}Gauss")
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_2)") 

#
# Gauss feps=1.e-6
#clear T Y err1 err2
Nphie = []
s=4
N=N0/s
# options(2) = nbre max iteration point fixe
options(2) = 15
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12 comme erreur absolue
options(3) = 1e-6
Ifail = []
for i=1:length(N),
    options(1) = N(i)
    [T,Y,nphie,ifail] = ode_gauss_v2(phi_vdp,[t0 tf],y0,options)
    Ifail = [Ifail length(find(ifail==-1))] 
    n=length(T)
    err1(i)=abs(y0(1)-Y(n,1))
    err2(i)=abs(y0(2)-Y(n,2))
    Nphie(i)=nphie
end
# 
figure(3) 
# vraie courbe
# ------------
subplot(1,2,1)
#hold on
Plots.plot(log10(Nphie),log10(err1))
i=round(length(N)/2)
#text(log10(Nphie(i)),log10(err1(i)),"\color{red}Gauss")
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_1)")     
subplot(1,2,2)
#hold on
Plots.plot(log10(Nphie),log10(err2))
i=round(length(N)/2)
#text(log10(Nphie(i)),log10(err2(i)),"\color{red}Gauss")
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_2)") 
#
# Gauss fpitermax=2
#clear T Y err1 err2
Nphie = []
s=4
N=N0/s
# options(2) = nbre max iteration point fixe
options(2) = 2
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12 comme erreur absolue
options(3) = 1e-12
Ifail = []
for i=1:length(N),
    options(1) = N(i)
    [T,Y,nphie,ifail] = ode_gauss_v2(phi_vdp,[t0 tf],y0,options)
    Ifail = [Ifail length(find(ifail==-1))] 
    n=length(T)
    err1(i)=abs(y0(1)-Y(n,1))
    err2(i)=abs(y0(2)-Y(n,2))
    Nphie(i)=nphie
end
figure(3) 
# vraie courbe
# ------------
subplot(1,2,1)
Plots.plot(log10(Nphie),log10(err1),"g")
i=round(length(N)/2)
#text(log10(Nphie(i)),log10(err1(i)),"\color{red}Gauss")
xlabel("log_{10}(nphie)")
ylabel("log_{10}(erreur pour y_1)")  
v=axis
axis([v(1) v(2) v(3)-2 v(4)])
subplot(1,2,2)
Plots.plot(log10(Nphie),log10(err2),"g")
i=round(length(N)/2)
#text(log10(Nphie(i)),log10(err2(i)),"\color{red}Gauss")
xlabel("log_{10}(fnphi)")
ylabel("log_{10}(erreur pour y_2)") 
v = axis
axis([v(1) v(2) v(3)-2 v(4)])
legend("gauss v2 feps=1.e-12","gauss v2 feps=1.e-6","gauss v2 fpitermax=2","Location","SouthWest")
#print -depsc fig_ordre_Gauss3
#
# Gauss + Newton
# ---------------
Nphie = []
s=2
N=N0/(2*s)
# options(2) = nbre max iteration point fixe
options(2) = 15
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12
# comme erreur absolue
options(3) = 1e-12
Ifail = []
for i=1:length(N),
    options(1) = N(i)
    [T,Y,nphie,ndphie,ifail] = ode_gauss_v3(phi_vdp,d_phi_vdp,[t0 tf],y0,options)
    Ifail = [Ifail length(find(ifail==-1))] 
    n=length(T)
    err1(i)=abs(y0(1)-Y(n,1))
    err2(i)=abs(y0(2)-Y(n,2))
    Nphie(i)=nphie+n*ndphie
end
# Inutile de trie car deja en ordre
#figure(1) # pour verifier l"ordre
# ---------------------
#subplot(1,2,1)
##hold on
#plot(log10(N0),log10(err1),"red")    
#subplot(1,2,2)
##hold on
#loglog(N0,err2,"red")
figure(2) 
# vraie courbe
# ------------
subplot(1,2,1)
Plots.plot(log10(Nphie),log10(err1),"red")
i=round(length(N)/2)
#text(log10(Nphie(i)),log10(err1(i)),"\color{red}Gauss")
#xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_1)")     
#v=axis
#axis([v(1) v(2) v(3)-4 v(4)])
subplot(1,2,2)
Plots.plot(log10(Nphie),log10(err2),"red")
i=round(length(N)/2)
#text(log10(Nphie(i)),log10(err2(i)),"\color{red}Gauss")
#xlabel("log_{10}(nphi)")
ylabel("log_{10}(erreur pour y_2)") 
#axis([v(1) v(2) v(3)-4 v(4)])
#legend("Euler", "Runge","Heun","RK4_1","RK4_2","gauss v2 feps=1.e-12","Location","SouthWest")