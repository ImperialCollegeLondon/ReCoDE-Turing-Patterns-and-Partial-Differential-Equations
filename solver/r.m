close all
clear all

x = readmatrix('IBVPx_2eqn_2D.dat');
y = readmatrix('IBVPy_2eqn_2D.dat');
Sol1 = readmatrix('IBVP1_2eqn_2D.dat');
Sol2 = readmatrix('IBVP2_2eqn_2D.dat');
nx = height(x(:,1))
ny = width(y(1,:))
nt = 2000
skip = 10

xr = max(x(:,1))
yr = max(y(1,:))

% Define function and time range
% Create animated plot


%syms u x1 y1
%u(x1,y1) = sin(x1*pi) + sin(y1*pi)

%Sol = Sol1 + Sol2    

for k = 1:skip:nt

    chr = int2str(k)

    f=figure(1);
    f.Position = [0 0 400 400]
    pcolor(x,y,Sol1(1+(k-1)*nx:nx+(k-1)*nx,1:ny))
    %pcolor(x,y,u(x,y))
    hold on
    colorbar
    title(append('Activator at time step: ',chr))
    shading interp
    caxis([0.95 1.05])
    colormap('cool')
    xlabel('x axis')
    ylabel('y axis')
    ylim([0 xr]);
    xlim([0 yr])
    exportgraphics(gcf,'testAnimated1.gif','Append',true);
    hold off
    
    f=figure(2);
    f.Position = [500 0 400 400]
    pcolor(x,y,Sol2(1+(k-1)*nx:nx+(k-1)*nx,1:ny))
    title(append('Substrate at time step: ',chr))
    %pcolor(x,y,u(x,y))
    hold on
    colorbar
    shading interp
    caxis([0.98 1.02])
    xlabel('x axis')
    ylabel('y axis')
    ylim([0 xr]);
    xlim([0 yr])
    exportgraphics(gcf,'testAnimated2.gif','Append',true);
    hold off



    %pause(0.05);
    clear('1')
    clear('2')
    hold off

    
end