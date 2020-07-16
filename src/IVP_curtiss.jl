using Plots
using LinearAlgebra

include("phi_curtiss.jl")
include("ode_euler_imp_v2.jl")

#
# ~gergaud/ENS/edo/Projet/ordre/IVP_curtiss.m
#
# Auteurs:  Joseph GERGAUD
# Date:     avril 2008
# Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
#           2; rue Camichel 31071 Toulouse FRANCE
# Email:    gergaud@enseeiht.fr
#***************************************************************************
#
# Intégration de l"equation differentiel de l"equation Curtiss et Hirschfelder
# ref: Hairer page 3 tome 2; pb raide
#
pyplot()
closeall()
y0=[0]
t0=0
tf=1.5
N0=[floor(tf*50/1.974) floor(tf*50/1.875) 50 100]; N=N0
option = zeros(3,1)
#=
for i=1:length(N)
  figure(i)
  #hold on
  # Euler
  T,Y=ode_euler(phi_toto,[t0 tf],y0,N[i])
  plot!(T,Y,LineWidth=2, label = "Euler")
  #
  # Runge
  T,Y=ode_runge(phi_curtiss,[t0 tf],y0,N[i])
  plot(T,Y,color:=black,LineWidth = 2,label = "Runge")

  #    
  # Heun
  T,Y=ode_heun(phi_curtiss,[t0 tf],y0,N[i],label = "Heun")
  plot(T,Y,color:=yellow,LineWidth = 2)

  # RK4 classique
  T,Y=ode_rk41(phi_curtiss,[t0 tf],y0,N[i],label = "RK41")
  plot(T,Y,color:=green,LineWidth=2)

  # RK4 regle 3/8
  [T,Y]=ode_rk42(phi_curtiss,[t0 tf],y0,N[i],label = "RK42")
  plot(T,Y,color:=blue,LineWidth = 2) 
  # Gauss
  option[2] = 40
  option[3] = 1e-6
  option[1] = N[i]
  T,Y,nphie,ifail=ode_gauss_v2(phi_curtiss,[t0 tf],y0,option)
  plot(T,Y,color:=red,LineWidth = 2,xlabel = "t",ylabel = "y(t)" , title = "Gauss");  
  
end
=#
#legend("Euler','Runge','Heun','RK41','RK42','Gauss")
#figure(5)
#hold on
#figure(6)
#hold on
plt = Plots.plot(layout = (2,2))
for i=1:length(N)
	#
	# Gauss
	
	option[2] = 40
	option[3] = 1e-6
	option[1] = N[i]
	#=
	T,Y,nphie,ifail=ode_gauss_v2(phi_curtiss,[t0 tf],y0,option)
	ifail
	subplot(2,2,i)
	Plots.plot!(plt[i],T,Y,color:=red,LineWidth=2);
	=#  
	T,Y,nphie,ifail=ode_euler_imp_v2(phi_curtiss,[t0 tf],y0,option)
	print(ifail)
	print("\n")
	#figure(6)
	#subplot(2,2,i)
	plot!(plt[i],T,Y,color=:green,LineWidth=2)
	#figure(5)
	#title("Gauss")
	#figure(6)
	#title("Euler implicite")
	# print("-depsc",fich)
end 
display(plt)
#=
for i=1:6
figure(i)
fich=["figcurtiss" int2str(i)]
print("-depsc",fich)
unix(["epstopdf figcurtiss' int2str(i) '.eps"])
end
=#
