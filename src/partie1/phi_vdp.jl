"""
 ~gergaud/ENS/Control/ODE/Matlab/euler.m

 Auteurs:  Joseph GERGAUD
 Date:     janvier 2000
 Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
           2, rue Camichel 31071 Toulouse FRANCE
 Email:    gergaud@enseeiht.fr
***************************************************************************

 deuxieme membre de l'equation differentiel de l'equation de Van der Pol
   ref: Hairer
"""
function phi_vdp(t, y)
    ypoint = [y[2]; (1 - y[1] * y[1]) * y[2] - y[1]]

    return ypoint
end