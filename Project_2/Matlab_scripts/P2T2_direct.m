clear;clc;close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model description                          %
% *******************************************%
% 7 DOF systems - 4 wheels, one rigid body   %
% pitching and rolling, pneumatic stiffness  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose : evaluating the response   %
% of a 7-DOF modelled vehicles to     %
% bumps on a flat road.               %
% Suspension and pneumatic stiffness  %
% are taken into consideration.       %
% The chassis is a simple rigid body. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Definition of inputs
% -1 assignment to enter the cicle. Commented value is the standard one.
ms = -1; % 1359;          %sprung masses [kg]
mn = -1; % 47;            %unsprung masses [kg]
Jx = -1; % 486;           %moment of inertia along x-axis [kg*m^2]
Jy = -1; % 2366;          %moment of inertia along y-axis [kg*m^2]
a = -1;  % 1.632/2;       %y distance between wheels and center of gravity [m]
b1 = -1; % 1.457;         %x distance between front wheels and front gravity [m]
b2 = -1; % 1.511;         %x distance between rear wheels and front gravity [m]
p = -1;  % 310000;        %tyre stiffness [N/m]
c = -1;  % 1450;          %suspension damping [Ns/m]
k0 = -1; % 30000;         %suspension stiffness [N/m]
k1 = -1; % 44000;         %anti-roll stiffness [N/m]
%*********************************************************************
% INSERT THE MAIN PARAMETERS RESPECTING A GIVEN RANGE
% If no value is inserted, the default one is used
%*********************************************************************
while ((ms<=800) || (ms>2000))
    fprintf('Leave blank to use the default value\n\n');
    ms=input('Insert the mass of the vehicle ( range 800-2000 [kg] ) = ');
    if isempty(ms)
        ms=1359;
    end
end
while ((mn<=15) || (mn>80))
    mn=input('Insert the mass of the wheels ( range 15-80 [kg] ) =');
    if isempty(mn)
        mn=47;
    end
end
while ((Jx<=200) || (Jx>1200))
    Jx=input('Insert the moment of inertia Jx ( range 200-1200 [kg*m^2] ) = ');
    if isempty(Jx)
        Jx=486;
    end
end
while ((Jy<=2000) || (Jy>12000))
    Jy=input('Insert the moment of inertia Jy ( range 2000-12000 [kg*m^2] ) = ');
    if isempty(Jy)
        Jy=2366;
    end
end
while ((a<=0.3) || (a>1.3))
    a=input('Insert y distance between wheels and center of gravity ( range 0.3-1.3 [m] ) = ');
    if isempty(a)
        a=1.632/2;
    end
end
while ((b1<=0.5) || (b1>2))
    b1=input('Insert x distance between front wheels and center of gravity ( range 0.5-2 [m] ) = ');
    if isempty(b1)
        b1=1.457;
    end
end
while ((b2<=0.5) || (b2>2))
    b2=input('Insert x distance between rear wheels and center of gravity ( range 0.5-2 [m] ) = ');
    if isempty(b2)
        b2=1.511;
    end
end
while ((p<=100000) || (p>400000))
    p=input('Insert the pneumatic stiffness ( range 100000-400000 [N/m] ) = ');
    if isempty(p)
        p=310000;
    end
end
while ((k0<=20000) || (k0>100000))
    k0=input('Insert the suspension stiffness ( range 20000-100000 [N/m] ) = ');
    if isempty(k0)
        k0=30000;
    end
end
while ((c<=1000) || (c>3000))
    c=input('Insert the suspension damping ( range 1000-3000 [Ns/m] ) = ');
    if isempty(c)
        c=1450;
    end
end
while ((k1<0) || (k1>600000))
    k1=input('Insert the anti-rollbar stiffness  ( max 600000 [kN/m], zero for no anti-roll ) = ');
    if isempty(k1)
        k1=44000;
    end
end

%% Definition of the system's matrices
M = diag([ms Jx Jy mn mn mn mn]);   %mass

% For ease of use, only the upper triangle of the matrix is defined, 
% and is then reflected upon the diagonal
C = [4*c  0      2*c*(b1*-b2)    -c      -c    -c     -c   ;
     0   4*a^2*c  0              -a*c     a*c   a*c   -a*c ;
     0    0      2*c*(b1^2+b2^2) -b1*c   -b1*c  b2*c   b2*c;
     0    0       0               c       0     0      0   ;
     0    0       0               0       c     0      0   ;
     0    0       0               0       0     c      0   ;
     0    0       0               0       0     0      c  ];    %damping
C = C + C' - diag(diag(C)); %make it symmetric
K = [4*k0 0              2*k0*(b1-b2)     -k0               -k0              -k0      -k0   ;
     0    k0*4*a^2       0                -a*k0              a*k0             a*k0    -a*k0 ;
     0    0              2*k0*(b1^2+b2^2) -b1*k0            -b1*k0            b2*k0    b2*k0;
     0    0              0                k0+p+k1           -k1               0        0    ;
     0    0              0                0                  k0+p+k1          0        0    ;
     0    0              0                0                  0                k0+p     0    ;
     0    0              0                0                  0                0        k0+p]; %stiffness
K = K + K' - diag(diag(K)); %make it symmetric


%% Definitions of the bumps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The use of mm in this section is not coherent with the rest of the code,%
% but is needed to have only integer numbers to use in the "for" loop.    %
% Inconcistencies are corrected afterwards.                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=(b1+b2)*1000;        %wheelbase [mm]
d=150;                 %height of the bump [mm]
L=l+2*d+d/3;           %distance between two bumps x=(1:1:L);

% Shaping of a single bump
y=zeros(1,L);          
for i=(l+1):(l+d)
    y(i)=y(i-1)+1;
end
for i=(l+d+1):(l+d+d/3)
    y(i)=d;
end
for i=(l+d+d/3+1):L
    y(i)=y(i-1)-1;
end

X=(0:1:10*L+L+L/2+50000)./1000;     %full length of the test, back to [m]

% Defining how each wheel "sees" all the bumps
Y1=[0,y,y,y,y,y,y,y,y,y,y,zeros(1,(L+L/2)),zeros(1,50000)];         %fr
Y2=[0,zeros(1,L/2),y,y,y,y,y,y,y,y,y,y,zeros(1,L),zeros(1,50000)];  %fl
Y3=[0,zeros(1,(L+L/2)),y,y,y,y,y,y,y,y,y,y,zeros(1,50000)];         %rl
Y4=[0,zeros(1,L),y,y,y,y,y,y,y,y,y,y,zeros(1,L/2),zeros(1,50000)];  %rr

acc=1;                              %acceleration [m/(s^2)]
T=sqrt((2*X)/acc);                  %time needed to complete the test

% Definition of some symbolic variables
syms t                                              %time
syms z(t) phi(t) theta(t) y1(t) y2(t) y3(t) y4(t)   %lagrangian coordinates
syms h1 h2 h3 h4                                    %bumps
syms sinbump(L,A,P)                                 %bumps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the use of ode45 for the direct resolution of the equations of      %
% motion, the trapezoidal bumps have been approximated with sinusoidal    %
% bumps (piecewise, symbolic). This is needed to minimize possible        %
% numerical errors and to improve speed, since Matlab's ode solvers do not%
% perform well with discontinuous functions.                              %
% The approximation introduces a small underestimation error.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sinbump(L,A,P) = A*(0.5+0.5*cos(2*pi*((t-P)/L)))*((P-L/2 < t) & (t < P+L/2));
% The "for" loops simply "scan" through the vector of the bumps and add a
% sine wave where needed.
h1=0;h2=0;h3=0;h4=0;
for ii=1:length(Y1)
    if Y1(ii) ~= 0 && Y1(ii-1) == 0
        st = T(ii); %save the time when a bump begins
        if h1 == 0
            first_bump = ii; %save the T index when the flat road ends, will be useful later
        end
    end
    if Y1(ii) ~= 0 && Y1(ii+1) == 0
        en = T(ii); %save the time when a bump ends
        h1 = h1 + sinbump(en-st,0.15,(en+st)/2); %add a sine bump
    end
end
for ii=1:length(Y2)
    if Y2(ii) ~= 0 && Y2(ii-1) == 0
        st = T(ii);
    end
    if Y2(ii) ~= 0 && Y2(ii+1) == 0
        en = T(ii);
        h2 = h2 + sinbump(en-st,0.15,(en+st)/2);
    end
end
for ii=1:length(Y3)
    if Y3(ii) ~= 0 && Y3(ii-1) == 0
        st = T(ii);
    end
    if Y3(ii) ~= 0 && Y3(ii+1) == 0
        en = T(ii);
        h3 = h3 + sinbump(en-st,0.15,(en+st)/2);
    end
end
for ii=1:length(Y4)
    if Y4(ii) ~= 0 && Y4(ii-1) == 0
        st = T(ii);
    end
    if Y4(ii) ~= 0 && Y4(ii+1) == 0
        en = T(ii);
        h4 = h4 + sinbump(en-st,0.15,(en+st)/2);
    end
end

%% Definition of the equations of motion
Y = [z; phi; theta; y1; y2; y3; y4]; %unknowns

F = [0; 0; 0; p*h1; p*h2; p*h3; p*h4]; %forces

odes = M*diff(Y,2) + C*diff(Y) + K*Y == F; %EoM

% Conversion of the symbolic equations into a form that ode45 can understand
[V,S] = odeToVectorField(odes); 
% V is the actual equation to solve, while S is a vector containing the
% substitutions made by odeToVectorField (which are not trivial). It will 
% be useful later, to know which components of the solution to plot.
FUN = matlabFunction(V,'vars', {'t','Y'});

%% Resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The test is divided into three parts:                                   %
% 1. Flat road before the bumps                                           %
% 2. Section with the bumps                                               %
% 3. Flat road after the bumps                                            %
% This is needed to minimize the numerical error, since Matlab's ode      % 
% solvers might have troubles when the solution changes suddenly.         % 
% This also allows to have an increased precision only where needed (the  % 
% bumps), thus improving speed.                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interval = [0 T(first_bump-1)];
yInit = zeros(1,14);
ySol_1 = ode45(FUN,interval,yInit);
opt = odeset('Refine',15,'RelTol',1e-10,'AbsTol',1e-10,'InitialStep',1e-3,'MaxStep',1e-2);
ySol_2 = odextend(ySol_1,FUN,T(end-50000),yInit,opt);
ySol_3 = odextend(ySol_2,FUN,T(end));

tValues_1 = linspace(0,T(first_bump-1),100);
tValues_2 = linspace(T(first_bump-1),T(end-50000),10000);
tValues_3 = linspace(T(end-50000),T(end),1000);

y1Values_1 = deval(ySol_1,tValues_1,1);
thetaValues_1 = deval(ySol_1,tValues_1,3);
y2Values_1 = deval(ySol_1,tValues_1,5);
y3Values_1 = deval(ySol_1,tValues_1,7);
y4Values_1 = deval(ySol_1,tValues_1,9);
zValues_1 = deval(ySol_1,tValues_1,11);
phiValues_1 = deval(ySol_1,tValues_1,13);
y1Values_2 = deval(ySol_2,tValues_2,1);
thetaValues_2 = deval(ySol_2,tValues_2,3);
y2Values_2 = deval(ySol_2,tValues_2,5);
y3Values_2 = deval(ySol_2,tValues_2,7);
y4Values_2 = deval(ySol_2,tValues_2,9);
zValues_2 = deval(ySol_2,tValues_2,11);
phiValues_2 = deval(ySol_2,tValues_2,13);
y1Values_3 = deval(ySol_3,tValues_3,1);
thetaValues_3 = deval(ySol_3,tValues_3,3);
y2Values_3 = deval(ySol_3,tValues_3,5);
y3Values_3 = deval(ySol_3,tValues_3,7);
y4Values_3 = deval(ySol_3,tValues_3,9);
zValues_3 = deval(ySol_3,tValues_3,11);
phiValues_3 = deval(ySol_3,tValues_3,13);

%% Plotting of the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
plot(X,Y1,'b-')
plot(X,Y2,'r-')
plot(X,Y3,'r--','LineWidth',1.5)
plot(X,Y4,'b--','LineWidth',1.5)
legend('front right (y1)','front left (y2)','rear left (y3)','rear right (y4)')
title('Profile of the bumps in space')
xlabel('Distance along the test strip [m]')
ylabel('Distance from the ground [mm]')
grid minor
set(gca, 'FontName', 'Helvetica')
set(gca, 'FontSize',15)
set(gca, 'GridAlpha',1)
grid on
set(gca, 'MinorGridAlpha',0.9)
set(gcf, 'Position', get(0, 'Screensize')-[0 0 10 10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
subplot(2,1,1)
hold on
plot(T,Y1,'b-')
plot(T,Y2,'r-')
plot(T,Y3,'r--','LineWidth',1.5)
plot(T,Y4,'b--','LineWidth',1.5)
legend('front right (y1)','front left (y2)','rear left (y3)','rear right (y4)')
title('Adjusted in time with a = 1 m/s^2')
xlabel('Time [s]')
ylabel('Distance from the ground [mm]')

grid minor
set(gca, 'FontName', 'Helvetica')
set(gca, 'FontSize',15)
set(gca, 'GridAlpha',1)
grid on
set(gca, 'MinorGridAlpha',0.9)

% since the solution has been split into three parts, a loop for legend 
% and graph representation is avoided not to mix the relative graph's colors
subplot(2,1,2)
hold on
plot(tValues_1,y1Values_1*1000,'b','Linewidth',1.5)
plot(tValues_1,y2Values_1*1000,'r','Linewidth',1.5)
plot(tValues_1,y3Values_1*1000,'g','Linewidth',1.5)
plot(tValues_1,y4Values_1*1000,'m','Linewidth',1.5)
plot(tValues_2,y1Values_2*1000,'b','Linewidth',1.5)
plot(tValues_2,y2Values_2*1000,'r','Linewidth',1.5)
plot(tValues_2,y3Values_2*1000,'g','Linewidth',1.5)
plot(tValues_2,y4Values_2*1000,'m','Linewidth',1.5)
plot(tValues_3,y1Values_3*1000,'b','Linewidth',1.5)
plot(tValues_3,y2Values_3*1000,'r','Linewidth',1.5)
plot(tValues_3,y3Values_3*1000,'g','Linewidth',1.5)
plot(tValues_3,y4Values_3*1000,'m','Linewidth',1.5)

% representation is in [mm] scale  of
% the four wheels' displacement is represented

xlabel('Time [s]')
legend('y1','y2','y3','y4')
ylabel('Displacement [mm]')
title('Wheel displacement')
grid minor
set(gca, 'FontName', 'Helvetica')
set(gca, 'FontSize',15)
set(gca, 'GridAlpha',1)
grid on
set(gca, 'MinorGridAlpha',0.9)
set(gcf, 'Position', get(0, 'Screensize')-[0 0 10 10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
hold on

% representation of the wheels' displacement restricted to the bumps

plot(tValues_2,y1Values_2*1000,'b','Linewidth',1.5)
plot(tValues_2,y2Values_2*1000,'r','Linewidth',1.5)

% representation is in [mm] scale 
% the four wheels' displacement is represented

xlabel('Time [s]')
legend('y1','y2')
ylabel('Displacement [mm]')
title('Wheel displacement rescaled, front wheels')
grid minor
set(gca, 'FontName', 'Helvetica')
set(gca, 'FontSize',15)
set(gca, 'GridAlpha',1)
grid on
set(gca, 'MinorGridAlpha',0.9)
set(gcf, 'Position', get(0, 'Screensize')-[0 0 10 10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
hold on

% representation of the wheels' displacement restricted to the bumps

plot(tValues_2,y3Values_2*1000,'m','Linewidth',1.5)
plot(tValues_2,y4Values_2*1000,'g','Linewidth',1.5)

% representation is in [mm] scale 
% the four wheels' displacement is represented

xlabel('Time [s]')
legend('y3','y4')
ylabel('Displacement [mm]')
title('Wheel displacement rescaled, rear wheels')
grid minor
set(gca, 'FontName', 'Helvetica')
set(gca, 'FontSize',15)
set(gca, 'GridAlpha',1)
grid on
set(gca, 'MinorGridAlpha',0.9)
set(gcf, 'Position', get(0, 'Screensize')-[0 0 10 10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
subplot(2,1,1);
hold on

% the rolling and pitch angle are represented in ° scale 

plot(tValues_1,phiValues_1*57.2958,'g','Linewidth',2)
plot(tValues_1,thetaValues_1*57.2958,'r','Linewidth',2)
plot(tValues_2,phiValues_2*57.2958,'g','Linewidth',2)
plot(tValues_2,thetaValues_2*57.2958,'r','Linewidth',2)
plot(tValues_3,phiValues_3*57.2958,'g','Linewidth',2)
plot(tValues_3,thetaValues_3*57.2958,'r','Linewidth',2)

legend({'\phi (Rolling)','\theta (Pitching)'});

xlabel('Time [s]')
ylabel('Angle [°]')
title('Rolling, Pitching')
grid minor
set(gca, 'FontName', 'Helvetica')
set(gca, 'FontSize',15)
set(gca, 'GridAlpha',1)
grid on
set(gca, 'MinorGridAlpha',0.9)
set(gcf, 'Position', get(0, 'Screensize')-[0 0 10 10]);

subplot(2,1,2);
time = [tValues_1 tValues_2 tValues_3];
theta = [thetaValues_1 thetaValues_2 thetaValues_3];
dtheta = diff(theta)./diff(time);
dtheta2 = diff(dtheta)./diff(time(1:end-1));
phi = [phiValues_1 phiValues_2 phiValues_3];
dphi = diff(phi)./diff(time);
dphi2 = diff(dphi)./diff(time(1:end-1));

hold on
plot(time(1:end-2),dphi2,'g','Linewidth',1.5)
plot(time(1:end-2),dtheta2,'r','Linewidth',1.5)
xlabel('Time [s]')
ylabel('Acceleration [rad/(s^2)]')
title('Body angular accelerations')
grid minor
set(gca, 'FontName', 'Helvetica')
set(gca, 'FontSize',15)
set(gca, 'GridAlpha',1)
grid on
set(gca, 'MinorGridAlpha',0.9)
set(gcf, 'Position', get(0, 'Screensize')-[0 0 10 10]);
legend({'\phi (Rolling)','\theta (Pitching)'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)
subplot(2,1,1);
hold on

% the position of the vehicle's center of gravity is represented in [mm] scale
plot(tValues_1,zValues_1*1000,'m','Linewidth',1.5)
plot(tValues_2,zValues_2*1000,'m','Linewidth',1.5)
plot(tValues_3,zValues_3*1000,'m','Linewidth',1.5)

xlabel('Time [s]')
ylabel('Displacement [mm]')
title('Center of gravity position')
grid minor
set(gca, 'FontName', 'Helvetica')
set(gca, 'FontSize',15)
set(gca, 'GridAlpha',1)
grid on
set(gca, 'MinorGridAlpha',0.9)
set(gcf, 'Position', get(0, 'Screensize')-[0 0 10 10]);
z = [zValues_1 zValues_2 zValues_3];
time = [tValues_1 tValues_2 tValues_3];
dz = diff(z)./diff(time);
dz2 = diff(dz)./diff(time(1:end-1));

subplot(2,1,2);
hold on
plot(time(1:end-2),dz2,'m','Linewidth',1.5)
xlabel('Time [s]')
ylabel('Acceleration [m/(s^2)]')
title('Center of gravity acceleration')
grid minor
set(gca, 'FontName', 'Helvetica')
set(gca, 'FontSize',15)
set(gca, 'GridAlpha',1)
grid on
set(gca, 'MinorGridAlpha',0.9)
set(gcf, 'Position', get(0, 'Screensize')-[0 0 10 10]);
