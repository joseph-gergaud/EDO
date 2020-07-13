"""
 ~gergaud/ENS/Control/ODE/Matlab/rk4_1.m

 Auteurs:  Joseph GERGAUD
 Date:     nov. 2005
 Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
           2, rue Camichel 31071 Toulouse FRANCE
 Email:    gergaud@enseeiht.fr
***************************************************************************

 fonction Phi pour "LE" schema de Runge-Kutta 4
"""
function rk4_1(f::Function, t, y, h) 
    k1 = f(t, y)
    k2 = f(t + h / 2, y + (h / 2) * k1)
    k3 = f(t + h / 2, y + (h / 2) * k2)
    k4 = f(t + h, y + h * k3)
    Phi=(1/6) * (k1 + 2*k2 + 2 * k3 + k4)

    return Phi
end
