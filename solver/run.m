close all
clear all

x = readmatrix('BVPx.dat');
y = readmatrix('BVPy.dat');
Sol1 = readmatrix('BVP1.dat');
Sol2 = readmatrix('BVP2.dat');


f=figure(1);
pcolor(x,y,Sol1)
hold on
shading interp
colorbar
hold off

f=figure(2);
pcolor(x,y,Sol2)
hold on
shading interp
colorbar
hold off