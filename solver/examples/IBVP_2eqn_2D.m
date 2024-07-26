%%%% This MATLAB script plots the solution to a coupled IBVP in two
%%%% dimensions

close all
clear all

% read in data

x = readmatrix('IBVPx_2eqn_2D.dat');
y = readmatrix('IBVPy_2eqn_2D.dat');
Sol1 = readmatrix('IBVP1_2eqn_2D.dat');
Sol2 = readmatrix('IBVP2_2eqn_2D.dat');
nx = height(x(:,1))
ny = width(y(1,:))
nt = 2000
skip = 20

xr = max(x(:,1))
yr = max(y(1,:))

% Create animated plot

for k = 1:skip:nt

    chr = int2str(k)

    f=figure(1);
    f.Position = [0 0 400 400]
    % be careful of indices
    pcolor(x,y,Sol1(1+(k-1)*nx:nx+(k-1)*nx,1:ny))
    hold on
    colorbar
    caxis([0.95 1.05])
    fontsize(f, 12, "points")
    title(append('Activator at time step: ',chr),'Interpreter','latex','FontSize',18)
    shading interp
   
    colormap('cool')
    xlabel('x axis','Interpreter','latex','FontSize',18)
    ylabel('y axis','Interpreter','latex','FontSize',18)

    ylim([0 xr]);
    xlim([0 yr])
    %exportgraphics(gcf,'examples/Activator.gif','Append',true);
    hold off
    
    f=figure(2);
    f.Position = [500 0 400 400]
    pcolor(x,y,Sol2(1+(k-1)*nx:nx+(k-1)*nx,1:ny))
    hold on
    colorbar
    caxis([0.98 1.02])
    colorbar
    shading interp
    fontsize(f, 12, "points")
    title(append('Substrate at time step: ',chr),'Interpreter','latex','FontSize',18)
    xlabel('x axis','Interpreter','latex','FontSize',18)
    ylabel('y axis','Interpreter','latex','FontSize',18)
    ylim([0 xr]);
    xlim([0 yr])
    %exportgraphics(gcf,'examples/Substrate.gif','Append',true);
    hold off

    pause(0.05);
    clear('1')
    clear('2')
    hold off
end