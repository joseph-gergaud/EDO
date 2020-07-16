#
# ~gergaud/ENS/edo/Projet/ordre/phi_curtiss.m
#
# Auteurs:  Joseph GERGAUD
# Date:     avril 2008
# Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
#           2; rue Camichel 31071 Toulouse FRANCE
# Email:    gergaud@enseeiht.fr
#***************************************************************************
#
# deuxieme membre de l"equation differentiel de l"equation de Curtiss et Hirschfelder
  # ref: Hairer page 2 tome 2
#
function phi_curtiss(t,y)
     ypoint=-50*(y.-cos(t))
     return ypoint
end
