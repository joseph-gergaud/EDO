
include("../src/ode_euler.jl")
include("../src/ode_runge.jl")
include("../src/ode_rk41.jl")
include("../src/ode_rk42.jl")
include("../src/ode_heun.jl")
include("../src/ode_euler_pasfixe.jl")

include("../src/ode_gauss_pf.jl")
include("../src/ode_euler_pf.jl")

include("../src/ode_gauss_newton.jl")

include("../src/d_fun_vdp.jl")
include("../src/fun_vdp.jl")

solutions = []
tol_erreur = 1e-3

@testset "test solvers" begin 

    t0 = 0
    tf = 2
    y0 = rand(2)
    phi(t,y) = [y[2]; (1 - y[1] * y[1]) * y[2] - y[1]]
    dphi(t,y) = [0 1; -2*y[1]*y[2]-1 1-y[1]^2]
    N = 10
    options = [N,15,1e-6]
    # euler
    T,Y = ode_euler(phi,[t0,tf],y0,N)
    @test isapprox(T, T, atol=tol_erreur)
    @test isapprox(Y, Y, atol=tol_erreur)
    # runge
    T,Y = ode_runge(phi,[t0,tf],y0,N)
    @test isapprox(T, T, atol=tol_erreur)
    @test isapprox(Y, Y, atol=tol_erreur)
    # rk41
    T,Y = ode_rk41(phi,[t0,tf],y0,N)
    @test isapprox(T, T, atol=tol_erreur)
    @test isapprox(Y, Y, atol=tol_erreur)
    # rk42
    T,Y = ode_rk42(phi,[t0,tf],y0,N)
    @test isapprox(T, T, atol=tol_erreur)
    @test isapprox(Y, Y, atol=tol_erreur)
    # heun
    T,Y = ode_heun(phi,[t0,tf],y0,N)
    @test isapprox(T, T, atol=tol_erreur)
    @test isapprox(Y, Y, atol=tol_erreur)
    # euler_pasfixe
    T,Y = ode_euler_pasfixe(phi,[t0,tf],y0,N)
    @test isapprox(T, T, atol=tol_erreur)
    @test isapprox(Y, Y, atol=tol_erreur)
    
    # gauss point fixe
    T,Y,nphie,ifail = ode_gauss_pf(phi,[t0,tf],y0,options)
    @test isapprox(T, T, atol=tol_erreur)
    @test isapprox(Y, Y, atol=tol_erreur)
    @test isapprox(nphie, nphie, atol=tol_erreur)
    @test isapprox(ifail, ifail, atol=tol_erreur)
    
    # euler point fixe
    T,Y,nphie,ifail = ode_euler_pf(phi,[t0,tf],y0,options)
    @test isapprox(T, T, atol=tol_erreur)
    @test isapprox(Y, Y, atol=tol_erreur)
    @test isapprox(nphie, nphie, atol=tol_erreur)
    @test isapprox(ifail, ifail, atol=tol_erreur)

    # gauss newton
    T,Y,nphie,ndphie,ifail = ode_gauss_newton(phi,dphi,[t0,tf],y0,options)
    @test isapprox(T, T, atol=tol_erreur)
    @test isapprox(Y, Y, atol=tol_erreur)
    @test isapprox(nphie, nphie, atol=tol_erreur)
    @test isapprox(ndphie, ndphie, atol=tol_erreur)
    @test isapprox(ifail, ifail, atol=tol_erreur)
end