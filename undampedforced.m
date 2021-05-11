xo = 2;
xdot = 0;
wn = 7;
w= 5;
Tn = 2*pi/wn;
Fo = 10;
t = 0:0.001:13;
beta = w/wn;
theta1 = wn.*t;
theta2 = w.*t;
k = 100;
zeta = 0.1;
wd = wn*sqrt(1-zeta^2);
x1 = xo.*cos(theta1);
x2 = (xdot/wn).*sin(theta1);
x3 = - (Fo/k).*(beta/(1-(beta)^2)).*sin(theta1);
x4 = (Fo/k).*(1/(1-beta^2)).*sin(theta2);
xundamped =x1+x2+x3+x4;

A = xo + (Fo/k)*((2*zeta*beta)/((1-beta^2)^2+(2*zeta*beta)^2));
B = (A*zeta*(wn/wd)) - (Fo/k)*(w/wd)*((1-beta^2)/((1-beta^2)^2+(2*zeta*beta)^2));
x5 = exp(-zeta*wn.*t);
x6 = A*cos(wd.*t) + B*sin(wd.*t);
x7 = x5.*x6;
x8 = (Fo/k)/((1-beta^2)^2+(2*zeta*beta)^2);
x9 = (1-beta^2).*sin(w.*t);
x10 = ((2*zeta*beta)^2).*cos(w.*t);
x11 = x8.*(x9 - x10);
xdamp = x7 + x11;

tiledlayout(2,1);
nexttile
plot(t,xundamped,'r','linewidth',2);
nexttile
plot(t,xdamp,'b','linewidth',2);
xlim([0,t(end)])
xlabel('time(t)');
ylabel('Response x(t)');
title(['Free spring-mass system(x(0) = ', num2str(xo),' m,  xdot(0) = ',num2str(xdot),' m/s,  wn = ',num2str(wn),' rad/s,  w = ',num2str(w), ' rad/s,  zeta = ',num2str(zeta)]);