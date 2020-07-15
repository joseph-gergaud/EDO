#
# ~gergaud/ENS/Control/ENS/edo/Projets/ordre/d_phi_vdp.m
#
# Auteurs:  Joseph GERGAUD
# Date:     december 2011
# Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
#           2; rue Camichel 31071 Toulouse FRANCE
# Email:    gergaud@enseeiht.fr
#***************************************************************************
#
# dérivée du deuxieme membre de l"equation differentiel de l"equation de Van der Pol
# ref: Hairer
#
function phi_dvp(t,y)
     dypoint=[0 1
	     -2*y[1]*y[2]-1 1-y[1]^2]

     return dypoint
end
