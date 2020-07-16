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
"""
y0=[2.00861986087484313650940188
     0]
t0=0
tf=6.6632868593231301896996820305
N0=10
N=N0
#
# Euler
for i=1:length(N)
     disp('Euler')
     [T,Y]=ode_euler(phi_vdp,[t0 tf],y0,N(i))
     [T Y]
     Y
     pause
          subplot(2,1,1)
          hold on
          plot(T,Y(:,1))
          subplot(2,1,2)
          hold on
          plot(T,Y(:,2))

     disp('Runge')
     [T,Y]=ode_runge(phi_vdp,[t0 tf],y0,N(i))
     [T Y]
     subplot(2,1,1)
     hold on
     plot(T,Y(:,1),'k')
     subplot(2,1,2)
     hold on
     plot(T,Y(:,2),'k')
     disp('Heun')
     [T,Y]=ode_heun(phi_vdp,[t0 tf],y0,N(i))
     [T Y]
     subplot(2,1,1)
     hold on
     plot(T,Y(:,1),'b')
     subplot(2,1,2)
     hold on
     plot(T,Y(:,2))
     disp('RK41')
     [T,Y]=ode_rk41(phi_vdp,[t0 tf],y0,N(i))
     [T Y]
     subplot(2,1,1)
     hold on
     plot(T,Y(:,1),'g')
     subplot(2,1,2)
     hold on
     plot(T,Y(:,2),'g')
     disp('RK42')
     [T,Y]=ode_rk42(phi_vdp,[t0 tf],y0,N(i))
     [T Y]
          subplot(2,1,1)
          hold on
          plot(T,Y(:,1),'r')
          subplot(2,1,2)
          hold on
          plot(T,Y(:,2),'r')
     # Gauss
     option(2) = 15
     # epsilon pour le test du point fixe
     # Il faut 1.e-12 car pour les valeurs de N grandes on atteint avec RK4 1.e-12 comme erreur absolue
     option(3) = 1e-6
     Ifail = []
     option(1) = N(i)
     option
     [T,Y,nphie,ifail]=ode_gauss_v2(phi_vdp,[t0 tf],y0,option)
     [T Y]
     nphie
     ifail
     Ifail = [Ifail length(find(ifail==-1))] 
end
