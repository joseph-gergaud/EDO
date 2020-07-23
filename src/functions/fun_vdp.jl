using Markdown
@doc doc"""
Deuxieme membre de l'equation differentiel de l'equation de Van der Pol

# Inputs 
  - t = l'instant t  
  - y = la valeur de y à cet instant

# Outputs
  - ypoint = ``\left\{ \begin{aligned}\begin{array}{c} \dot{y}_{1}(t)= &y_{2}(t) \\ \dot{y}_{2}(t)= &(1-y_{1}^{2}(t)) y_{2}(t)-y_{1}(t) \end{array} \end{aligned}\right.`` 
""" 
function fun_vdp(t, y)
    ypoint = [y[2]; (1 - y[1] * y[1]) * y[2] - y[1]]
    return ypoint
end