using Markdown
@doc doc"""
dérivée du deuxieme membre de l'equation differentiel de l'equation de Van der Pol
# Inputs 
  - t = l'instant t  
  - y = la valeur de y à cet instant

# Outputs
  - dypoint = ``\left[\begin{array}{cc} 0 & 1 \\ -2 y_1(t)y_2(t)-1 & 1-y_1(t)^2 \end{array} \right]``
"""
function d_fun_vdp(t,y)
     
    dypoint = [0 1; -2*y[1]*y[2]-1 1-y[1]^2]

    return dypoint
end
