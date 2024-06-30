close all
clear all

X1 = readmatrix('IVBP_1.dat');
X2 = readmatrix('IVBP_2.dat');
w = width(X1(1,:))
h = height(X1(:,1))

x = X1(2:h,1)
t = X1(1,2:w)
sol = X1(2:h,2:w)

sol2 = X2(2:h,2:w)

syms u x1 t1
u(x1,t1) = exp(-4*pi^2*t1)*sin(2*pi*x1)

% Define function and time range
% Create animated plot
f=figure(1);
    

for k = 1:15:length(t)
    plot(x, sol(:, k),'LineWidth',2,'color','k');
    hold on
    plot(x, sol(:, 1),'LineWidth',2,'color','k','LineStyle','--');
    plot(x, sol2(:, 1),'LineWidth',2,'color','c','LineStyle','--');
    
    plot(x, sol2(:, k),'LineWidth',2,'color','c','LineStyle','-');
    %plot(x, u(x,t(k)),'LineWidth',1,'color','c');
    xlabel('x');
    ylabel('y');
    %title(sprintf('t = %0.1f', t(k)));
    ylim([0 2]);
    pause(0.08);
    clear('1')
    hold off
end