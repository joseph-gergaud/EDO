"""
 ~gergaud/ENS/edo/Projets/ordre/ode_gauss_3.m

 Auteurs:  Joseph GERGAUD
 Date:     december 2011
 Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
           2, rue Camichel 31071 Toulouse FRANCE
 Email:    gergaud@enseeiht.fr
***************************************************************************

 Numerical integration by Gauss method of order 4
 The solution on the non linear equation is computed by Newton method

 function [T,Y,nphie,ndphie,ifail,KK]=ode_gauss_v3(phi,dphi,t0tf,y0,option)
 Input parameters
 ----------------
 phi = second member ypoint=phi(t,y)
 t0tf = [t0,tf]
 y0 = initial point
 option(1) = N = number of step
 option(2) = fpitermax = maximum number of iterations for the fixed point
 option(3) = fpeps = epsilon for the test of progress in the fixed point

 Output parameters
 -----------------
 T = vector of time
 Y = Matrix of solution
 The line i of [T Y] contains ti and y(ti)
 ifail(i) = 1 = computation successful for the fixed point on [t_i,t_{i+1}] 
 ifail(i) = -1 = computation failed for the fixed point on [t_i,t_{i+1}]: maximum number of iteration
 is attained in the fixed point
 nphie = number of evaluation of phi

 Local variables
 ---------------
 c1, c2, a11, a12, a21, a22, b1 and b2 = coefficient of Butcher array
 N = number of steps
 h = constant step
 k1 and k2 = k1 and k2 of the Runge-Kutta scheme
 delta1y
"""
function ode_gauss_v3(phi::Function,dphi::Function,t0tf,y0,option)

    # Initialisation
    # --------------
    # Coefficent of the Butcher's array
    c1 = 1/2 - sqrt(3)/6
    c2 = 1/2 + sqrt(3)/6
    a11 = 1/4 
    a12 = 1/4 - sqrt(3)/6
    a21 = 1/4 + sqrt(3)/6
    a22 = 1/4
    b1 = 1/2 
    b2 = 1/2
    I = [1 0;0 1]
    
    # Input variables
    t = t0tf[1]
    tf = t0tf[2]
    N = Int(option[1])
    fpitermax = option[2]
    fpeps = option[3]

    # Output variables
    ifail = ones(N) 
    nphie = 0 
    ndphie = 0
    T = zeros(N+1)
    T[1] = t
    y = y0[:]
    n = length(y0)
    Y = zeros(N+1,n)
    Y[1,:] = y'
    # Local variables
    h = (tf - t)/N
    KK = []
    
    # Boucle principale
    for i=1:N,
        t1 = t + c1*h
        t2 = t + c2*h
        K = [phi(t1,y) ;phi(t2,y)]
        if KK == []
            KK = K
        else
            KK = [KK K]
        end
        k1 = K[1:n]
        k2 = K[n+1:2*n]
        nphie = nphie + 2
        normprog = 1 
        nbiter = 0
        while (normprog > fpeps) & (nbiter <= fpitermax)     # Newton
            g1 = y + h * (a11*k1 + a12*k2)
            g2 = y + h * (a21*k1 + a22*k2)
            phik1 = phi(t1,g1)
            phik2 = phi(t2,g2)
            FK0 = [k1-phik1  k2-phik2]
            nphie = nphie + 2
            d_phik1 = dphi(t1,g1)
            d_phik2 = dphi(t2,g2)
            ndphie = ndphie + 2
            JFK0 = [I-a11*h*d_phik1 -a12*h*d_phik1; -a21*h*d_phik2  I-a22*h*d_phik2]
            println(size(k1))
            println(size(FK0))
            deltaK = -JFK0 \ FK0
            normprog = norm(deltaK)
            K = K + deltaK
            KK = [KK K]
            k1 = K[1:n] 
            k2 = K[n+1:2*n]
            nbiter = nbiter + 1        
        end # end while

        if (normprog > fpeps)
            ifail[i] = -1
        else
            ifail[i] = nbiter
        end
        
        t = t+h
        T[i+1] = t
        y = y + h*(b1*k1 + b2*k2)
        Y[i+1,:] = y'	
    end # end for
    T = T[:]
    
    return T,Y,nphie,ndphie,ifail,KK
end

