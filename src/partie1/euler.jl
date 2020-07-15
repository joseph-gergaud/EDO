#
# ~gergaud/ENS/Control/ODE/Matlab/euler.m
#
# Auteurs:  Joseph GERGAUD
# Date:     nov. 2005
# Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
#           2; rue Camichel 31071 Toulouse FRANCE
# Email:    gergaud@enseeiht.fr
#***************************************************************************
#
# fonction Phi pour le schema d'Euler.
#
function euler(f::Function,t,y,h)
  Phi=f(t,y,h)

  return Phi
end
