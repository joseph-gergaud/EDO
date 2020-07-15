#
# ~gergaud/ENS/Control/ODE/Matlab/heun.m
#
# Auteurs:  Joseph GERGAUD
# Date:     nov. 2005
# Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
#           2; rue Camichel 31071 Toulouse FRANCE
# Email:    gergaud@enseeiht.fr
#***************************************************************************
#
# fonction Phi pour le schema de Runge
#
function heun(f::Function,t,y,h)
  k1=f(t,y)
  k2=f(t+h/3,y+(h/3)*k1)
  k3=f(t+(2*h/3),y+(2*h/3)*k2)
  Phi=(1/4)*(k1+3*k3)
  return Phi
end
