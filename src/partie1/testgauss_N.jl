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

y0=[2.008619860874843136, 0]
n = length(y0)
t0=0
tf=6.663286859323130189
N0=[120:60:1080 1200:600:10800] # multiple de 12 pour avoir des nombres entier si on divise par 4 ou 3
N=N0
#
# Gauss Newton

#clear T Y err1 err2
Nphie = []
s=2
N=N0/(2*s)
N=10
# option(2) = nbre max iteration point fixe
option(2) = 15
# epsilon pour le test du point fixe
# Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12
# comme erreur absolue
option(3) = 1e-12
Ifail = []
for i=1:length(N),
option(1) = N(i)
[T,Y,nphie,ndphie,ifail]=ode_gauss_v3(phi_vdp,d_phi_vdp,[t0 tf],y0,option)
Ifail = [Ifail length(find(ifail==-1))] 
n=length(T)
err1(i)=abs(y0(1)-Y(n,1))
err2(i)=abs(y0(2)-Y(n,2))
Nphie(i)=nphie+n*ndphie
end
# Inutile de trie car deja en ordre
figure(1) # pour verifier l"ordre
# ---------------------
subplot(1,2,1)
#hold on
plot(log10(N0),log10(err1),"r")    
subplot(1,2,2)
#hold on
loglog(N0,err2,"r")
figure(2) 
# vraie courbe
# ------------
subplot(1,2,1)
plot(log10(Nphie),log10(err1),"r")
i=round(length(N)/2)
#text(log10(Nphie(i)),log10(err1(i)),"\color{red}Gauss")
xlabel("log_{10}(fe)")
ylabel("log_{10}(erreur pour y_1)")     
v=axis
axis([v(1) v(2) v(3)-4 v(4)])
subplot(1,2,2)
plot(log10(Nphie),log10(err2),"r")
i=round(length(N)/2)
#text(log10(Nphie(i)),log10(err2(i)),"\color{red}Gauss")
xlabel("log_{10}(nphi)")
ylabel("log_{10}(erreur pour y_2)") 
axis([v(1) v(2) v(3)-4 v(4)])
#legend("Euler", "Runge","Heun","RK4_1","RK4_2","gauss v2 feps=1.e-12","Location","SouthWest")
figure(3) 


[T,Y,nphie,ndphie,ifail,KK] = ode_gauss_v3(phi_vdp,d_phi_vdp,[t0 tf],y0,option)
