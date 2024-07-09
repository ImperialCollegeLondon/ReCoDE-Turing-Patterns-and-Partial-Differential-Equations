close all
clear all

BVP = readmatrix('BVP.dat');
%y = readmatrix('BVPy.dat');
%Sol1 = readmatrix('BVP1.dat');
%Sol2 = readmatrix('BVP2.dat');


epsi = 0.0001
epsi2 = 1/sqrt(epsi)

x = BVP(1:height(BVP(:,1)),1)
y = BVP(1:height(BVP(:,1)),3)

syms g x1
g(x1)= exp(2/sqrt(epsi))*x1 + 2*exp((1-x1)/sqrt(epsi)) - 2*exp((x1+1)/sqrt(epsi)) -x1
g(x1) = g(x1)/(1-exp(2/sqrt(epsi)))

f=figure(1);

plot(x,y,'LineStyle','-','LineWidth',3)
hold on
fplot(g,[0 1],'LineStyle','--','LineWidth',3,'color','r')
fontsize(f, 22, "points")
t=title('$\epsilon = 0.0001$','Interpreter','latex')
ylabel('u')
xlabel('x')
legend('numerical solution','analytical solution','location','northwest','fontsize',20)
hold off
