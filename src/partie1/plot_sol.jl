using Plots
function plot_sol(T,Y,c,labels)
    # plot les solutions
    # T = temps
    # Y = solutions
    # c = couleur
    plot!(T,Y[:,1], color=c, xlabel="t", ylabel="y_1(t)",label=labels[1], subplot=1)    
    plot!(T,Y[:,2], color=c, xlabel="t", ylabel="y_2(t)",label=labels[2], subplot=2)    
    plot!(Y[:,1],Y[:,2], color=c, xlabel="y_1(t)", ylabel="y_2(t)",label=labels[3], subplot=3)
    display(plt)
end