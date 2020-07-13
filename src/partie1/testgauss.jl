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
clear all; close all;
unix('rm testgauss.txt')
diary('testgauss.txt')
y0=[2.008619860874843136
    0];
t0=0;
tf=6.663286859323130189;
N0=[120:60:1080 1200:600:10800]; # multiple de 12 pour avoir des nombres entier si on divise par 4 ou 3
N=N0;
N=10
#
# Gauss point fixe
N= 10;
option(2) = 15;
option(3) = 1e-6;
option(1) = N;
[T,Y,nphie,ifail,KK] = ode_gauss_v2(phi_vdp,[t0 tf],y0,option);
#
disp('Point fixe')
disp('----------')
disp('[T Y]')
[T Y]
nphie
ifail
I = find(ifail==-1);
ifail(I) = option(2) ;

# KK(:,cumsum(ifail+1))
#
# Gauss Newton
[T,Y,nphie,ndphie,ifail,KK] = ode_gauss_v3(phi_vdp,d_phi_vdp,[t0 tf],y0,option);
disp('Newton')
disp('------')
disp('[T Y]')
[T Y]
nphie
ndphie
ifail
KK

subplot(2,2,1)
plot(T,Y(:,1))
xlabel('t')
ylabel('y_1(t)')
subplot(2,2,2)
plot(T,Y(:,2))
xlabel('t')
ylabel('y_2(t)')
subplot(2,2,3)
plot(Y(:,1),Y(:,2))
ylabel('y_1(t)')
ylabel('y_2(t)')

N= 200;
option(2) = 15;
option(3) = 1e-6;
option(1) = N;
#[T,Y,nphie,ifail] = ode_gauss_v2(phi_vdp,[t0 tf],y0,option);
#
subplot(2,2,1)
hold on
plot(T,Y(:,1),'r')
xlabel('t')
ylabel('y_1(t)')
subplot(2,2,2)
hold on
plot(T,Y(:,2),'r')
xlabel('t')
ylabel('y_2(t)')
subplot(2,2,3)
hold on
plot(Y(:,1),Y(:,2),'r')
xlabel('y_1(t)')
ylabel('y_2(t)')

diary
