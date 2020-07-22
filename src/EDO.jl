module EDO

using LaTeXStrings
using Test
using LinearAlgebra
using Markdown

include("ode_gauss_newton.jl")  ; export ode_gauss_newton
include("ode_gauss_pf.jl")      ; export ode_gauss_pf
include("ode_gauss.jl")         ; export ode_gauss
include("ode_euler.jl")         ; export ode_euler
include("ode_euler_pasfixe.jl") ; export ode_euler_pasfixe
include("ode_runge.jl")         ; export ode_runge
include("ode_rk41.jl")          ; export ode_rk41
include("ode_rk42.jl")          ; export ode_rk42
include("ode_heun.jl")          ; export ode_heun
include("plot_sol.jl")          ; export plot_sol
include("d_fun_vdp.jl")         ; export d_fun_vdp
include("fun_vdp.jl")           ; export fun_vdp

end