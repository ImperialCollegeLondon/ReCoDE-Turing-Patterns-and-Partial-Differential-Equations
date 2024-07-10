close all
clear all

X1 = readmatrix('IVBP_1eqn_1D.dat');
%X2 = readmatrix('IVBP2_2eqn_1D.dat');
w = width(X1(1,:))
h = height(X1(:,1))

x = X1(2:h,1)
t = X1(1,2:w)
sol = X1(2:h,2:w)
%sol2 = X2(2:h,2:w)

%sol2 = X2(2:h,2:w)

syms u x1 t1
u(x1,t1) = exp(-4*pi^2*t1)*sin(2*pi*x1)

% Define function and time range
% Create animated plot
f=figure(1);

skip = 1

for k = 1:skip:70
    chr = int2str(k)
    plot(x, sol(:, k),'LineWidth',2,'color','k');
    hold on
    plot(x, sol(:, 1),'LineWidth',2,'color','k','LineStyle','--');

    fontsize(f, 12, "points")
    title(append('Solution at time step: ',chr),'Interpreter','latex','FontSize',18)

    %plot(x, sol2(:, 1),'LineWidth',2,'color','c','LineStyle','--');
    
    %plot(x, sol2(:, k),'LineWidth',2,'color','c','LineStyle','-');
    %plot(x, u(x,t(k)),'LineWidth',1,'color','c');
    xlabel('x axis','Interpreter','latex','FontSize',18)
    ylabel('u axis','Interpreter','latex','FontSize',18)
    %title(sprintf('t = %0.1f', t(k)));
    %ylim([0.98 1.02]);
        exportgraphics(gcf,'examples/non_linear.gif','Append',true);
    %pause(0.2)
    clear('1')
    hold off
end

skip = 150

for k = 71:skip:length(t)
    chr = int2str(k)
    plot(x, sol(:, k),'LineWidth',2,'color','k');
    hold on
    plot(x, sol(:, 1),'LineWidth',2,'color','k','LineStyle','--');

    fontsize(f, 12, "points")
    title(append('Solution at time step: ',chr),'Interpreter','latex','FontSize',18)

    %plot(x, sol2(:, 1),'LineWidth',2,'color','c','LineStyle','--');
    
    %plot(x, sol2(:, k),'LineWidth',2,'color','c','LineStyle','-');
    %plot(x, u(x,t(k)),'LineWidth',1,'color','c');
    xlabel('x axis','Interpreter','latex','FontSize',18)
    ylabel('u axis','Interpreter','latex','FontSize',18)
    %title(sprintf('t = %0.1f', t(k)));
    %ylim([0.98 1.02]);
        exportgraphics(gcf,'examples/non_linear.gif','Append',true);
    %pause(0.2)
    clear('1')
    hold off
end