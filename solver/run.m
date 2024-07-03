close all
clear all

x = readmatrix('BVPx.dat');
y = readmatrix('BVPy.dat');
Sol1 = readmatrix('BVP1.dat');
Sol2 = readmatrix('BVP2.dat');


syms u x1 y1
u(x1,y1) = y1*(1-y1)*x1^3
C = (x.^2).*(y.^2)
%C = y*(1-y).*x.^3 - (x.^5)/10 + x/10



f=figure(1);
surf(x,y,Sol1)
hold on
shading interp
colorbar
ylabel('y')
hold off



f=figure(2);
surf(x,y,C)
hold on

shading interp
colorbar
hold off