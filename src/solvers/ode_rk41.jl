"""
Heun method 

# Description
    Numerical integration of the CauchX[i-1,:]'s problem
    x_point(t) = phi(t,x(t))
    x(t_0) = x_0

# Usage
    T, X = ode_rk4(phi,t0tf,X[i-1,:]0,N)

# Inputs
    f    - function     : second member of the ode whith the interface 
                          xpoint = phi(t, x)
                              t    - real     : time,
                              x = vector of R^n with the same dimension of x0
    t0tf - real(2)      : intial and final time  [t0,tf]
    x0   - real(n)      : initial point
    N    - integer      : number of steps (>1)
    
# Outputs
    T    - real(N+1,1)  : vector of times
    X    - real(N+1,n)  : Matrix of solution, The line i of [T Y] contains ti and x_i
"""
function ode_rk41(phi::Function,t0tf,x0,N)

    N = Int(N)
    n = length(x0)
    T = zeros(N+1,1)
    X = zeros(N+1,n)
    
    h = (t0tf[2]-t0tf[1])/N
    T[1] = 0
    X[1,:] = x0  
    
    for i = 2:N+1
        k1 = phi(T[i-1], X[i-1,:])
        k2 = phi(T[i-1] + h / 2, X[i-1,:] + (h / 2) * k1)
        k3 = phi(T[i-1] + h / 2, X[i-1,:] + (h / 2) * k2)
        k4 = phi(T[i-1] + h, X[i-1,:] + h * k3)
        X[i,:] = X[i-1,:] + (h/6) * (k1 + 2*k2 + 2 * k3 + k4)
    end
    return T,X
end   