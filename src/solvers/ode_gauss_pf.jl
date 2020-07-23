using LinearAlgebra
"""
Numerical integration by Gauss method of order 4

# Usage 
    - [T,Y,nphie,ifail,KK]=ode_gauss_pf(phi,t0tf,y0,options)

# Input parameters
    - phi = second member ypoint=phi(t,y)
    - t0tf = [t0,tf]
    - y0 = initial point
    - options(1) = N = number of step
    - options(2) = fpitermax = maximum number of iterations for the fixed point
    - options(3) = fpeps = epsilon for the test of progress in the fixed point

# Output parameters
    - T = vector of time
    - Y = Matrix of solution, The line i of [T Y] contains ti and y(ti)
    - ifail(i) = 1 = computation successful for the fixed point on [t_i,t_{i+1}] 
    - ifail(i) = -1 = computation failed for the fixed point on [t_i,t_{i+1}]: maximum number of iteration is attained in the fixed point
    - nphie = number of evaluation of phi

"""
function ode_gauss_pf(phi::Function, t0tf, y0, options)

    # Initialisation
    # --------------
    # Coefficents of the Butcher's array
    c1 = 1/2 - sqrt(3)/6  ; c2 = 1/2 + sqrt(3)/6
    a11 = 1/4             ; a12 = 1/4 - sqrt(3)/6
    a21 = 1/4 + sqrt(3)/6 ; a22 = 1/4
    b1 = 1/2              ; b2 = 1/2
    
    # Input variables
    t = t0tf[1];  tf = t0tf[2]
    N = Int(options[1])
    fpitermax = options[2]
    fpeps = options[3]

    # Output variables
    ifail = ones(N) 
    nphie = 0
    T = zeros(N + 1)
    T[1] = t
    y = y0[:]
    n = length(y0)
    Y = zeros(N + 1,n)
    Y[1,:] = y'

    # Local variables
    h = (tf - t)/N
       
    # Boucle principale
    for i=1:N
        t1 = t + c1 * h
        t2 = t + c2 * h
        k1 = phi(t1, y)
        k2 = phi(t2, y)
        nphie = nphie + 2
        normprog = 1 
        nbiter = 1
        while (normprog > fpeps) & (nbiter <= fpitermax) # point fixe
            k1new = phi(t1,y+h*(a11*k1 + a12*k2))
            k2new = phi(t2,y+h*(a21*k1 + a22*k2))
            normprog = norm([k1new-k1 k2new-k2])
            k1 = k1new
            k2 = k2new
            nbiter = nbiter + 1
            nphie = nphie + 2            
        end 
        if (normprog > fpeps)
            ifail[i] = -1
        else
            ifail[i] = nbiter - 1
        end
        t = t + h
        T[i + 1] = t
        y = y + h*(b1*k1 + b2*k2)
        Y[i + 1,:] = y'	
    end 
    return T, Y, nphie, ifail
end