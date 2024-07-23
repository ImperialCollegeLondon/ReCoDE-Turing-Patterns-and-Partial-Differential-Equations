%%%% This MATLAB script plots the solution to a coupled IBVP in one
%%%% dimension


close all
clear all

% read in a set data

X1 = readmatrix('IVBP1_2eqn_1D.dat');
X2 = readmatrix('IVBP2_2eqn_1D.dat');
w = width(X1(1,:))
h = height(X1(:,1))

x = X1(2:h,1)
t = X1(1,2:w)
sol1 = X1(2:h,2:w)
sol2 = X2(2:h,2:w)

% plot - animation

f=figure(1);

skip = 1

for k = 1:skip:length(t)
    chr = int2str(k)
    plot(x, sol1(:, k),'LineWidth',2,'color','b');
    hold on
    plot(x, sol2(:, k),'LineWidth',2,'color','r','LineStyle','-');
    plot(x, sol1(:, 1),'LineWidth',2,'color','b','LineStyle','--');
    fontsize(f, 12, "points")
    title(append('Predator/Prey Densities at time step: ',chr),'Interpreter','latex','FontSize',18)
    plot(x, sol2(:, 1),'LineWidth',2,'color','r','LineStyle','--');
    legend('Prey','Predator','location','northwest','fontsize',20,'Interpreter','latex')
    xlabel('$x$ axis','Interpreter','latex','FontSize',18)
    ylabel('$u$ and $v$ axis','Interpreter','latex','FontSize',18)
    %exportgraphics(gcf,'examples/predator_prey.gif','Append',true);
    hold off
    pause(0.2)
    clear('1')
end
