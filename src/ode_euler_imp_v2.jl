"""
 ~gergaud/ENS/edo/Matlab/ordre_un_pas/ode_euler_imp_v2.m

 Auteurs:  Joseph GERGAUD
 Date:     mars 2008
 Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
           2; rue Camichel 31071 Toulouse FRANCE
 Email:    gergaud@enseeiht.fr
***************************************************************************

 Numerical integration by implicit Euler

 function [T,Y,nphie,ifail]=ode_euler_imp_v2(phi,t0tf,y0,option)
 Input parameters
 ----------------
 phi = second member ypoint=phi[t,y]
 t0tf = [t0,tf]
 y0 = initial point
 option[1] = N = number of step
 option[2] = fpitermax = maximum number of iterations for the fixed point
 option[3] = fpeps = epsilon for the test of progress in the fixed point

 Output parameters
 -----------------
 T = vector of time
 Y = Matrix of solution
 The line i of [T Y] contains ti & y[ti]
 ifail[i] = 1 = computation successful for the fixed point on [t_i,t_[i+1]] 
 ifail[i] = -1 = computation failed for the fixed point on [t_i,t_[i+1]]: maximum number of iteration
 is attained in the fixed point
 nphie = number of evaluation of phi

 Local variables
 ---------------
 c1; c2; a11; a12; a21; a22; b1 & b2 = coefficient of Butcher array
 N = number of steps
 h = constant step
 k1 and k2 = k1 & k2 of the Runge-Kutta scheme
 delta1y
"""
function ode_euler_imp_v2(phi::Function,t0tf,y0,option)

    # Initialisation
    # --------------
    # Coefficent of the Butcher's array
    c = 1
    a = 1
    b = 1
    #
    # Input variables
    t = t0tf[1]; tf = t0tf[2]
    N = option[1]
    N = Int(N)
    T = zeros(N+1,1)
    Y = zeros(N+1,length(y0))
    fpitermax = option[2]
    fpeps = option[3]
    # Output variables
    ifail = ones(N,1); nphie = 1
    T[1] = t
    y = y0[:]
    Y[1,:] = y'
    # Local variables
    h = (tf - t)/N
    #
    # Boucle principale
    for i=1:N
    t1 = t + c*h
    k1 = phi(t1,y)
    normprog = 1; nbiter = 0
    while (normprog > fpeps) && (nbiter < fpitermax)     # fixed point
        k1new = phi(t1,y+h*a*k1)
        normprog = norm(k1new-k1)
        k1 = k1new
        nbiter = nbiter + 1
        nphie = nphie + 1
    end                                                  # end while
    if (normprog > fpeps)
        ifail[i] = -1
    else
        ifail[i] = 1
    end
    t = t+h
    T[i+1] = t
    y = y + h*(b*k1)
    Y[i+1,:] = y';	
    end                           # end for
    T = T[:]
    return T,Y,nphie,ifail

end
