function plot_sol(T,Y,c);
    % plot les soltions
    % T = temps;
    % Y = solutions
    % c = couleur
    subplot(2,2,1)
    hold on
    plot(T,Y(:,1),c)
    xlabel('t')
    ylabel('y_1(t)')
    subplot(2,2,2)
    hold on
    plot(T,Y(:,2),c)
    xlabel('t')
    ylabel('y_2(t)')
    subplot(2,2,3:4)
    hold on
    plot(Y(:,1),Y(:,2),c)
    xlabel('y_1(t)')
    ylabel('y_2(t)')
return;