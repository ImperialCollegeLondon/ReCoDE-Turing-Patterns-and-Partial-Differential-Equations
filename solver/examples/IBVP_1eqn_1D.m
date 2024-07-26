%%%% This MATLAB script plots the solution to a single IBVP in one
%%%% dimension

close all
clear all

% read in data
X1 = readmatrix('IVBP_1eqn_1D.dat');
w = width(X1(1,:))
h = height(X1(:,1))
x = X1(2:h,1)
t = X1(1,2:w)
sol = X1(2:h,2:w)


% analytical solution
syms u x1 t1
u(x1,t1) = exp(-4*pi^2*t1)*sin(2*pi*x1)

f=figure(1);

skip = 1

% animation of plot
for k = 1:skip:length(t)
    chr = int2str(k)
    plot(x, sol(:, k),'LineWidth',2,'color','b');
    hold on
    plot(x, sol(:, 1),'LineWidth',2,'color','b','LineStyle','--');
    fontsize(f, 12, "points")   
    legend('Prey','Predator','location','northeast','fontsize',20)
    xlabel('x axis','Interpreter','latex','FontSize',18)
    ylabel('u axis','Interpreter','latex','FontSize',18)
    hold off
    pause(0.2)
    clear('1')
end

