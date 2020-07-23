module EDO

using LaTeXStrings
using Test
using LinearAlgebra
using Markdown
using JLD2

include("solvers/ode_gauss_newton.jl")  ; export ode_gauss_newton
include("solvers/ode_gauss_pf.jl")      ; export ode_gauss_pf
include("solvers/ode_gauss.jl")         ; export ode_gauss
include("solvers/ode_euler.jl")         ; export ode_euler
include("solvers/ode_euler_pasfixe.jl") ; export ode_euler_pasfixe
include("solvers/ode_runge.jl")         ; export ode_runge
include("solvers/ode_rk41.jl")          ; export ode_rk41
include("solvers/ode_rk42.jl")          ; export ode_rk42
include("solvers/ode_heun.jl")          ; export ode_heun
include("plot_sol.jl")                  ; export plot_sol
include("functions/d_fun_vdp.jl")       ; export d_fun_vdp
include("functions/fun_vdp.jl")         ; export fun_vdp

end