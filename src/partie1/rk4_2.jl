"""
 ~gergaud/ENS/Control/ODE/Matlab/rk4_2m

 Auteurs:  Joseph GERGAUD
 Date:     nov. 2005
 Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
           2, rue Camichel 31071 Toulouse FRANCE
 Email:    gergaud@enseeiht.fr
***************************************************************************

 fonction Phi pour le schema de Runge-Kutta 4 regle des 3/8
"""
function rk4_2(f::Function, t, y, h) 
    
    k1 = f(t, y)
    k2 = f(t + h / 3, y + (h / 3) * k1)
    k3 = f(t + 2 * h / 3, y + h*(-(1/3) * k1 + k2))
    k4 = f(t + h, y + h * (k1 - k2 + k3))

    Phi = (1 / 8) * (k1 + 3 * k2+ 3 * k3 + k4)

    return Phi
end