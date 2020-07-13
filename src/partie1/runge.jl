"""
 ~gergaud/ENS/Control/ODE/Matlab/runge.m

 Auteurs:  Joseph GERGAUD
 Date:     nov. 2005
 Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
           2, rue Camichel 31071 Toulouse FRANCE
 Email:    gergaud@enseeiht.fr
***************************************************************************

 fonction Phi pour le schema de Runge
"""
function runge(f::Function, t, y, h)
    
    k1 = y + (h / 2) * f(t, y)
    Phi = f(t + h / 2, k1)

    return Phi
end