%%%% This MATLAB script plots the solution to a single BVP in one dimension

close all
clear all

% Read in data and set x and y
BVP = readmatrix('BVP.dat');
x = BVP(1:height(BVP(:,1)),1);
y = BVP(1:height(BVP(:,1)),3);

% Symbolic solution to compare
syms g x1
epsi = 0.01;
g(x1)= 1 + x1 + (3-1-1)*(exp(x1/epsi)-1)/(exp(1/epsi)-1);


f=figure(1);
plot(x,y,'LineStyle','-','LineWidth',3)
hold on
fplot(g,[0 1],'LineStyle','--','LineWidth',3,'color','r')
fontsize(f, 22, "points")
t=title('','Interpreter','latex')
ylabel('u','Interpreter','latex')
xlabel('x','Interpreter','latex')
legend('numerical solution','analytical solution','location','northwest','fontsize',20)
hold off
