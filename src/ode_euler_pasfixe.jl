"""
~gergaud/ENS/edo/Projet/ordre/ode_euler.m
Auteurs:  Joseph GERGAUD 

Date:     avril 2008

Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
          2, rue Camichel 31071 Toulouse FRANCE
Email:    gergaud@enseeiht.fr

**************************************************************************
Programme d'integration par le schema d'Euler a pas fixe et constant

"""
function ode_euler_pasfixe(f::Function,t0tf,y0,N)
    t0 = t0tf[1]
    tf = t0tf[2]
    y = y0
    N = Int(N)
    T = zeros(N + 1)
    Y = zeros(N + 1,length(y0))
    T[1] = t0
    Y[1,:] = y0'
    h = (tf-t0)/N
    t = t0     
    for i=1:N
        t=t+h
        T[i+1]=t
        y = y + h * f(t,y)
        Y[i+1,:] = y'
    end
    T = T' 
    return T,Y
end