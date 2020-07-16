"""
 ~gergaud/ENS/edo/Projet/ordre/IVP_curtiss.m

 Auteurs:  Joseph GERGAUD
 Date:     december 2011
 Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
           2, rue Camichel 31071 Toulouse FRANCE
 Email:    gergaudenseeiht.fr
***************************************************************************

  pb raide
"""

include("phi_stiff.jl")
y0 = [10]
t0 = 0
tf = 1.5
#N0=[fix(tf*50/1.974) fix(tf*50/1.875) 50 100] 
#N=N0
N = [30 40 80 100]#

for i=1:length(N)
    figure(i)
    hold on
    # Euler
    [T,Y] = ode_euler(phi_stiff,[t0 tf],y0,N(i))
    plot(T,Y,'LineWidth',2)
    hold on
    #
    # Runge
    [T,Y] = ode_runge(phi_stiff,[t0 tf],y0,N(i))
    plot(T,Y,'k','LineWidth',2)

    #    
    # Heun
    [T,Y] = ode_heun(phi_stiff,[t0 tf],y0,N(i))
    plot(T,Y,'y','LineWidth',2)

    # RK4 classique
    [T,Y] = ode_rk41(phi_stiff,[t0 tf],y0,N(i))
    plot(T,Y,'g','LineWidth',2)

    # RK4 regle 3/8
    [T,Y] = ode_rk42(phi_stiff,[t0 tf],y0,N(i))
    plot(T,Y,'m','LineWidth',2) 
    # Gauss
    option(2) = 40
    option(3) = 1e-6
    option(1) = N(i)
    [T,Y,nphie,ifail] = ode_gauss_v2(phi_stiff,[t0 tf],y0,option)
    plot(T,Y,'r','LineWidth',2)  
    ifail
    xlabel('t')
    ylabel('y(t)')
end
    legend('Euler','Runge','Heun','RK41','RK42','Gauss')
    #figure(5)
    #hold on
    #figure(6)
    #hold on
    #for i=1:length(N),
    #
    # Gauss
    #option(2) = 40
    #option(3) = 1e-6
    #option(1) = N(i)
    #[T,Y,nphie,ifail] = ode_gauss_v2(phi_stiff,[t0 tf],y0,option)
    #ifail
    #figure(5)
    #subplot(2,2,i)
    #plot(T,Y,'r','LineWidth',2)  
    #[T,Y,nphie,ifail] = ode_euler_imp_v2(phi_stiff,[t0 tf],y0,option)
    #ifail
    #figure(6)
    #subplot(2,2,i)
    #plot(T,Y,'g','LineWidth',2)
    #figure(5)
    #title('Gauss')
    #figure(6)
    #title('Euler implicite')
    # print('-depsc',fich)
    #end 

for i=1:4
    figure(i)
    fich=['figstiff' int2str(i)]
    print('-depsc',fich)
    unix(['epstopdf figstiff' int2str(i) '.eps'])
end
