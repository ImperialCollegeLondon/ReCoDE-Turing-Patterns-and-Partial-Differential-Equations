close all
clear all

X = readmatrix('IVBP_2.dat');
w = width(X(1,:))
h = height(X(:,1))

x = X(2:h,1)
t = X(1,2:w)
sol = X(2:h,2:w)

syms u x1 t1
u(x1,t1) = exp(-4*pi^2*t1)*sin(2*pi*x1)

% Define function and time range
% Create animated plot
f=figure(1);
for k = 1:10:length(t)
    plot(x, sol(:, k),'LineWidth',2,'color','k');
    hold on
    %plot(x, u(x,t(k)),'LineWidth',1,'color','c');
    xlabel('x');
    ylabel('y');
    %title(sprintf('t = %0.1f', t(k)));
    %ylim([-1.5, 1.5]);
    pause(0.1);
    %clear('1')
    %hold off
end