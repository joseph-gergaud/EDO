"""
Deuxieme membre de l"equation différentielle de l'équation de Curtiss et Hirschfelder
ref: Hairer page 2 tome 2
"""
function fun_curtiss(t,y)
    ypoint = -50*(y.-cos(t))
    return ypoint
end
