using Plots

"""
 ~gergaud/ENS/Control/ODE/Matlab/test_gauss.m

 Auteurs:  Joseph GERGAUD
 Date:     nov. 2005
 Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
           2, rue Camichel 31071 Toulouse FRANCE
 Email:    gergaudenseeiht.fr
***************************************************************************

 Int�gration de l'equation differentiel de l'equation de Van der Pol
 ref: Hairer
"""

include("ode_gauss_v2.jl")
include("ode_gauss_v3.jl")
include("phi_vdp.jl")
include("d_phi_vdp.jl")
#rm("testgauss.txt")
#diary("testgauss.txt")
y0=[2.008619860874843136; 0]
t0=0
tf=6.663286859323130189
N0=[120:60:1080; 1200:600:10800] # multiple de 12 pour avoir des nombres entier si on divise par 4 ou 3
N=N0
N=10
#
# Gauss point fixe
N= 10
options = zeros(3)
options[1] = N
options[2] = 15
options[3] = 1e-6
T,Y,nphie,ifail,KK = ode_gauss_v2(phi_vdp,[t0 tf],y0,options)
#
println("Point fixe")
println("----------")
println("[T Y]")
println([T Y])
println(nphie)
println(ifail)
I = findall(ifail==-1)
ifail[I] .= options[2]

# KK(:,cumsum(ifail+1))
#
# Gauss Newton
T,Y,nphie,ndphie,ifail,KK = ode_gauss_v3(phi_vdp,d_phi_vdp,[t0 tf],y0,options)
println("Newton")
println("------")
println("[T Y]")
println([T Y])
println(nphie)
println(ndphie)
println(ifail)
println(KK)

pyplot()
plt = Plots.plot(layout=(2,2))

Plots.plot!(T,Y[:,1],xlabel="t",ylabel="y_1(t)",subplot=1)
Plots.plot!(T,Y[:,2],xlabel="t",ylabel="y_2(t)",subplot=2)
Plots.plot!(Y[:,1],Y[:,2],xlabel="y_1(t)",ylabel="y_2(t)",subplot=3)

N= 200
options[2] = 15
options[3] = 1e-6
options[1] = N

[T,Y,nphie,ifail] = ode_gauss_v2(phi_vdp,[t0 tf],y0,options)

Plots.plot!(T,Y[:,1],xlabel="t",ylabel="y_1(t)",subplot=1)
Plots.plot!(T,Y[:,2],xlabel="t",ylabel="y_2(t)",subplot=2)
Plots.plot!(Y[:,1],Y[:,2],xlabel="y_1(t)",ylabel="y_2(t)",subplot=3)
