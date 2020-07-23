
include("../src/solvers/ode_euler.jl")
include("../src/solvers/ode_runge.jl")
include("../src/solvers/ode_rk41.jl")
include("../src/solvers/ode_rk42.jl")
include("../src/solvers/ode_heun.jl")
include("../src/solvers/ode_euler_pasfixe.jl")

include("../src/solvers/ode_gauss_pf.jl")
include("../src/solvers/ode_euler_pf.jl")

include("../src/solvers/ode_gauss_newton.jl")

include("../src/functions/d_fun_vdp.jl")
include("../src/functions/fun_vdp.jl")

solutions = []
tol_erreur = 1e-3

@testset "test avec l'eq de VDP" begin 

    # le premier cas de test 
    y0 = [2.008619860874843136, 0]
    t0 = 0
    tf = 6.663286859323130189
    phi(t,y)  = [y[2]; (1 - y[1] * y[1]) * y[2] - y[1]]
    dphi(t,y) = [0 1; -2*y[1]*y[2]-1 1-y[1]^2]
    N = 25
    options = [N,15,1e-6]
    @load "test_vdp.jld2"

    
    # euler
    T,Y = ode_euler(phi,[t0,tf],y0,N)
    #T_euler,Y_euler = ode_euler(phi,[t0,tf],y0,N)
    @test isapprox(T, T_euler, atol=tol_erreur)
    @test isapprox(Y, Y_euler, atol=tol_erreur)
    
    # runge
    T,Y = ode_runge(phi,[t0,tf],y0,N)
    #T_runge,Y_runge = ode_runge(phi,[t0,tf],y0,N)
    @test isapprox(T, T_runge, atol=tol_erreur)
    @test isapprox(Y, Y_runge, atol=tol_erreur)
    
    # rk41
    T,Y = ode_rk41(phi,[t0,tf],y0,N)
    #T_rk41,Y_rk41 = ode_rk41(phi,[t0,tf],y0,N)
    @test isapprox(T, T_rk41, atol=tol_erreur)
    @test isapprox(Y, Y_rk41, atol=tol_erreur)
    
    # rk42
    T,Y = ode_rk42(phi,[t0,tf],y0,N)
    #T_rk42,Y_rk42 = ode_rk42(phi,[t0,tf],y0,N)
    @test isapprox(T, T_rk42, atol=tol_erreur)
    @test isapprox(Y, Y_rk42, atol=tol_erreur)
    
    # heun
    T,Y = ode_heun(phi,[t0,tf],y0,N)
    #T_heun,Y_heun = ode_heun(phi,[t0,tf],y0,N)
    @test isapprox(T, T_heun, atol=tol_erreur)
    @test isapprox(Y, Y_heun, atol=tol_erreur)
    
    # euler_pasfixe
    T,Y = ode_euler_pasfixe(phi,[t0,tf],y0,N)
    #T_euler_pasfixe,Y_euler_pasfixe = ode_euler_pasfixe(phi,[t0,tf],y0,N)
    @test isapprox(T, T_euler_pasfixe, atol=tol_erreur)
    @test isapprox(Y, Y_euler_pasfixe, atol=tol_erreur)
    
    # gauss point fixe
    T,Y,nphie,ifail = ode_gauss_pf(phi,[t0,tf],y0,options)
    #T_gauss_pf,Y_gauss_pf,nphie_gauss_pf,ifail_gauss_pf = ode_gauss_pf(phi,[t0,tf],y0,options)
    @test isapprox(T, T_gauss_pf, atol=tol_erreur)
    @test isapprox(Y, Y_gauss_pf, atol=tol_erreur)
    @test isapprox(nphie, nphie_gauss_pf, atol=tol_erreur)
    @test isapprox(ifail, ifail_gauss_pf, atol=tol_erreur)
    
    # euler point fixe
    T,Y,nphie,ifail = ode_euler_pf(phi,[t0,tf],y0,options)
    #T_euler_pf,Y_euler_pf,nphie_euler_pf,ifail_euler_pf = ode_euler_pf(phi,[t0,tf],y0,options)
    @test isapprox(T, T_euler_pf, atol=tol_erreur)
    @test isapprox(Y, Y_euler_pf, atol=tol_erreur)
    @test isapprox(nphie, nphie_euler_pf, atol=tol_erreur)
    @test isapprox(ifail, ifail_euler_pf, atol=tol_erreur)

    # gauss newton
    T,Y,nphie,ndphie,ifail = ode_gauss_newton(phi,dphi,[t0,tf],y0,options)
    #T_newton,Y_newton,nphie_newton,ndphie_newton,ifail_newton = ode_gauss_newton(phi,dphi,[t0,tf],y0,options)
    @test isapprox(T, T_newton, atol=tol_erreur)
    @test isapprox(Y, Y_newton, atol=tol_erreur)
    @test isapprox(nphie, nphie_newton, atol=tol_erreur)
    @test isapprox(ndphie, ndphie_newton, atol=tol_erreur)
    @test isapprox(ifail, ifail_newton, atol=tol_erreur)

    #@save "test_vdp.jld2" T_euler Y_euler T_runge Y_runge T_rk41 Y_rk41 T_rk42 Y_rk42 T_heun Y_heun T_euler_pasfixe Y_euler_pasfixe T_gauss_pf Y_gauss_pf nphie_gauss_pf ifail_gauss_pf T_euler_pf Y_euler_pf nphie_euler_pf ifail_euler_pf T_newton Y_newton nphie_newton ndphie_newton ifail_newton
end