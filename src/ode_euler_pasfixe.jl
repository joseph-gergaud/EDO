"""
Programme d'integration par le schema d'Euler à pas constant

# Usage
    - T,Y,nphie,ifail=ode_euler_pasfixe(phi,t0tf,y0,N)

# Input parameters 
    - phi = second member ypoint=phi[t,y]
    - t0tf = [t0,tf]
    - y0 = initial point
    - N = number of step

# Output parameters
    - T = vector of time
    - Y = Matrix of solution, The line i of [T Y] contains ti & y[ti]
"""
function ode_euler_pasfixe(phi::Function,t0tf,y0,N)
        
    t0 = t0tf[1]; tf = t0tf[2]    
    y = y0      ; t = t0
    N = Int(N)  ; T = zeros(N + 1)
    Y = zeros(N + 1,length(y0))
    T[1] = t0   ; Y[1,:] = y0'
    h = (tf-t0)/N         
    for i=1:N
        t=t+h
        T[i+1]=t
        y = y + h * phi(t,y)
        Y[i+1,:] = y'
    end
    T = T' 
    return T,Y
end