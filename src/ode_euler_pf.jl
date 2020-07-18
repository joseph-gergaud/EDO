"""
Numerical integration by implicit Euler

# Usage
    - T,Y,nphie,ifail=ode_euler_pf(phi,t0tf,y0,options)

# Input parameters 
    - phi = second member ypoint=phi[t,y]
    - t0tf = [t0,tf]
    - y0 = initial point
    - options[1] = N = number of step
    - options[2] = fpitermax = maximum number of iterations for the fixed point
    - options[3] = fpeps = epsilon for the test of progress in the fixed point

# Output parameters
    - T = vector of time
    - Y = Matrix of solution, The line i of [T Y] contains ti & y[ti]
    - ifail[i] = 1 = computation successful for the fixed point on [t_i,t_[i+1]] 
    - ifail[i] = -1 = computation failed for the fixed point on [t_i,t_[i+1]]: maximum number of iteration is attained in the fixed point
    - nphie = number of evaluation of phi
"""
function ode_euler_pf(phi::Function,t0tf,y0,options)

    # Initialisation
    # --------------
    # Coefficent of the Butcher's array
    c = 1
    a = 1
    b = 1
    
    # Input variables
    t = t0tf[1]; tf = t0tf[2]
    N = options[1]
    N = Int(N)
    T = zeros(N+1,1)
    Y = zeros(N+1,length(y0))
    fpitermax = options[2]
    fpeps = options[3]

    # Output variables
    ifail = ones(N,1)
    nphie = 0
    T[1] = t
    y = y0[:]
    Y[1,:] = y'

    # Local variables
    h = (tf - t)/N
    
    # Boucle principale
    for i=1:N
        t1 = t + c*h
        k1 = phi(t1,y)
        nphie = nphie + 1
        normprog = 1
        nbiter = 0
        while (normprog > fpeps) && (nbiter < fpitermax)     # point fixe
            k1new = phi(t1,y+h*a*k1)
            nphie = nphie + 1
            normprog = norm(k1new-k1)
            k1 = k1new
            nbiter = nbiter + 1
        end
        if (normprog > fpeps)
            ifail[i] = -1
        else
            ifail[i] = 1
        end
        t = t+h
        T[i+1] = t
        y = y + h*(b*k1)
        Y[i+1,:] = y'	
    end
    return T,Y,nphie,ifail
end