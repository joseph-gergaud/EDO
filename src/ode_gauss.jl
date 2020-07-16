using LinearAlgebra
"""
 ~gergaud/ENS/edo/Matlab/ordre_un_pas/odegauss.m

 Auteurs:  Joseph GERGAUD
 Date:     mars 2008
 Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
           2; rue Camichel 31071 Toulouse FRANCE
 Email:    gergaud@enseeiht.fr
***************************************************************************

 Numerical integration by Gauss method of order 4

 function [T,Y,nphie,ifail]=ode_gauss(phi,t0tf,y0,option)
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
function ode_gauss(phi::Function,t0tf,y0,option)

    # Initialisation
    # --------------
    # Coefficent of the Butcher's array
    c1 = 1/2 - sqrt(3)/6  ; c2 = 1/2 + sqrt(3)/6
    a11 = 1/4             ; a12 = 1/4 - sqrt(3)/6
    a21 = 1/4 + sqrt(3)/6 ; a22 = 1/4
    b1 = 1/2              ; b2 = 1/2
    #
    # Input variables
    t = t0tf[1]; tf = t0tf[2]
    N = option[1]
    fpitermax = option[2]
    fpeps = option[3]
    # Output variables
    ifail = ones(N)
    k1 = []
    k2 = []
    nphie = 1
    T = zeros(N + 1)
    T[1] = t
    y = y0[:]
    T = zeros(N + 1,length(y))
    Y[1,:] = y'
    # Local variables
    h = (tf - t)/N
    #
    # Boucle principale
    for i=1:N
        t1 = t + c1*h
        t2 = t + c2*h
        k = phi(t,y)
        delta1y = h*c1*k            # c1=a11+a12
        delta2y = h*c2*k            # c2=a21+a22
        normprog = 1; nbiter = 0
        while (normprog > fpeps) && (nbiter <= fpitermax)     # fixed point
            k1 = phi(t1,y+delta1y)
            k2 = phi(t2,y+delta2y)
            newdelta1y = h*(a11*k1 + a12*k2)
            newdelta2y = h*(a21*k1 + a22*k2)
            normprog = norm([newdelta1y-delta1y  newdelta2y-delta2y]')
            delta1y = newdelta1y
            delta2y = newdelta2y
            nbiter = nbiter + 1
            nphie = nphie + 2
        end                                                  # end while
        if (normprog > fpeps)
            ifail[i] = -1
        else
            ifail[i] = 1
        end
        t = t+h
        T[i+1] = t
        y = y + h*(b1*k1 + b2*k2)
        Y[i+1,:] = y'	
    end                           # end for
    T = T[:]
    return T,Y,nphie,ifail
end
