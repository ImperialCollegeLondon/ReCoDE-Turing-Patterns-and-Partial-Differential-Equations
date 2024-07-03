close all
clear all

x = readmatrix('IBVPx.dat');
y = readmatrix('IBVPy.dat');
Sol1 = readmatrix('IBVP1.dat');
Sol2 = readmatrix('IBVP2.dat');
nx = height(x(:,1))
ny = width(y(1,:))
nt = 2000
skip = 1

xr = max(x(:,1))
yr = max(y(1,:))

% Define function and time range
% Create animated plot


syms u x1 y1
u(x1,y1) = sin(x1*pi) + sin(y1*pi)

    
for k = 1:skip:nt
    f=figure(1);
    f.Position = [0 0 500 500]
    pcolor(x,y,Sol1(1+(k-1)*nx:nx+(k-1)*nx,1:ny))
    %pcolor(x,y,u(x,y))
    hold on
    colorbar
    shading interp
    caxis([0 1])
    ylim([0 xr]);
    xlim([0 yr])
    pause(0.2);
    clear('1')
    hold off

    
end