#=========================================================================================
#
#    Heun method 
#
#    Description
#
#        Numerical integration of the Cauchy's problem
#        x_point(t) = f(t,x(t))
#        x(t_0) = x_0
#
#-------------------------------------------------------------------------------------------
#
#    Usage
#
#        T, X = ode_rk4(f,t0tf,y0,N)
#
#    Inputs
#        f    - function     : second member of the ode whith the interface 
#                              xpoint = f(t, x)
#                                  t    - real     : time,
#                                  x = vector of R^n with the same dimension of x0
#        t0tf - real(2)      : intial and final time  [t0,tf]
#        x0   - real(n)      : initial point
#        N    - integer      : number of steps (>1)
#
#    Outputs
#        T    - real(N+1,1)  : vector of times
#        X    - real(N+1,n)  : Matrix of solution
#        The line i of [T Y] contains ti and x_i
#
=###############################################################################################
include("rk4_1.jl")
include("rk4_2.jl")
function ode_rk4(f::Function,t0tf,x0,N,phi::Int)

    N = Int(N)
    n = length(x0)
    T = zeros(N+1,1)
    X = zeros(N+1,n)
    
    h = (t0tf[2]-t0tf[1])/N
    T[1] = 0
    X[1,:] = x0  
    
    for i = 2:N+1         
        #T[i] = (i-1)*h
        #k1 = f(T[i-1],X[i-1,:])
        #k2 = f(T[i-1] + h / 2, X[i-1,:] + h*(1/2)*k1)
        #k3 = f(T[i-1] + h / 2, X[i-1,:] + h*(0*k1 + (1/2)*k2))
        #k4 = f(T[i-1] + h,X[i-1,:] + h*(0*k1 + 0*k2 + k3))
        #X[i,:] = X[i-1,:] + h*((1/6)*k1 +(2/6)*k2 + (2/6)*k3 + (1/6)*k4)
        if phi == 2
            X[i,:] = X[i-1,:] + h*rk4_2(f,T[i-1],X[i-1,:],h)
        else
            X[i,:] = X[i-1,:] + h*rk4_1(f,T[i-1],X[i-1,:],h)
        end
    end

    return T,X
end
    
    




    