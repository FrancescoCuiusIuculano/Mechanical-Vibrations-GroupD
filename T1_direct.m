clear;clc;close all

ms = 1839;
mn = 40;
Jx = 606;
Jy = 4192;
a = 1.515/2;
b1 = 1.3;
b2 = 1.4;
p = 250;
c = 2500;
k0 = 10000;
k1 = 600;
syms t
syms z(t) phi(t) theta(t) y1(t) y2(t) y3(t) y4(t)
syms h1(t) h2(t) h3(t) h4(t)

h1=sin(t);h2=cos(t);h3=cos(t+1);h4=sin(t+1);
ode1 = ms*diff(z,2) ...
     + c*(4*diff(z)-diff(y1)-diff(y2)-diff(y3)-diff(y4)+2*(-b1+b2)*diff(theta)) ...
     + k0*(4*z-y1-y2-y3-y4+2*(-b1+b2)*theta) == 0;
ode2 = Jx*diff(phi,2) ...
     + c*a*(-diff(y1)+diff(y2)+diff(y3)-diff(y4)+4*a*diff(phi)) ...
     + k0*a*(-y1+y2+y3-y4+4*a*phi) + k1*(phi-(y1-y2)/(2*a)) == 0;
ode3 = Jy*diff(theta,2) ...
     + c*(-b1*(2*diff(z)-diff(y1)-diff(y2)-2*b1*diff(theta)) ...
          +b2*(2*diff(z)-diff(y3)-diff(y4)+2*b2*diff(theta))) ...
     + k0*(-b1*(2*z-y1-y2-2*b1*theta) + b2*(2*z-y3-y4+2*b2*theta)) == 0;
ode4 = mn*diff(y1,2) ...
     - c*(diff(z)-diff(y1)+a*diff(phi)-b1*diff(theta)) ...
     - k0*(z-y1+a*phi-b1*theta) ...
     - (k1/(2*a))*(phi-(y1-y2)/(2*a)) + p*(y1-h1) == 0;
ode5 = mn*diff(y2,2) ...
     - c*(diff(z)-diff(y2)-a*diff(phi)-b1*diff(theta)) ...
     - k0*(z-y2-a*phi-b1*theta) ...
     - (k1/(2*a))*(phi-(y1-y2)/(2*a)) + p*(y2-h2) == 0;
ode6 = mn*diff(y3,2) ...
     - c*(diff(z)-diff(y3)-a*diff(phi)+b2*diff(theta)) ...
     - k0*(z-y3-a*phi+b2*theta) + p*(y3-h3) == 0;
ode7 = mn*diff(y4,2) ...
     - c*(diff(z)-diff(y4)+a*diff(phi)+b2*diff(theta)) ...
     - k0*(z-y4+a*phi+b2*theta) + p*(y4-h4) == 0;
odes = [ode1 ode2 ode3 ode4 ode5 ode6 ode7];

[V,S] = odeToVectorField(odes);
M = matlabFunction(V,'vars', {'t','Y'});

interval = [0 200];
yInit = zeros(1,14);
ySol = ode45(M,interval,yInit);
tValues = linspace(0,200,10000);

y1Values = deval(ySol,tValues,1);
thetaValues = deval(ySol,tValues,3);
y2Values = deval(ySol,tValues,5);
y3Values = deval(ySol,tValues,7);
y4Values = deval(ySol,tValues,9);
zValues = deval(ySol,tValues,11);
phiValues = deval(ySol,tValues,13);

hold on
plot(tValues,y1Values)
plot(tValues,y2Values)
plot(tValues,y3Values)
plot(tValues,y4Values)
plot(tValues,zValues)
plot(tValues,phiValues)
plot(tValues,thetaValues)

legend('y1','y2','y3','y4','z','\phi','\theta')