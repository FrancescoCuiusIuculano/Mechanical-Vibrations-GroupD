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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Since T doesn't have a constant step, it is "translated" into a linspace%
% t, which has the same length, but constant step.                        %
% The vectors representing the bumps are then re-calculated accordingly.  %
% This passage is needed in oreder to get the correct result with the     %
% convolution.                                                            %
% Even though this approach might seem a bit more twisted than needed, it %
% has been chosen in order to allow easier and faster manipulation of the %
% input data.                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = linspace(0,T(end),100000);
h1 = interp1(T,Y1,t);
h2 = interp1(T,Y2,t);
h3 = interp1(T,Y3,t);
h4 = interp1(T,Y4,t);

%% Matrices in the state space 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A*diff(Y) + B*Y == F, where Y=[diff(q);q]                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [zeros(7,7), M;
        M, C];
    
B = [-M, zeros(7,7);
        zeros(7,7), K];

F = [zeros(10,length(t)); p*h1; p*h2; p*h3; p*h4]; %forces

[V, D] = eig(A\B); %calculation of eigenvectors and eigenvalues

Mt = transpose(V)*A*V;

Kt = transpose(V)*B*V;

% Normalizing the modal matrix
Vn= zeros(14);            %preallocating for optimizing speed
for j =1:14
    nom = 1/sqrt(Mt(j,j));
    for i = 1:14
        Vn(i,j)=V(i,j)*nom;
    end
end

Ft = transpose(Vn)*F; % Force in the modal space

%% Resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The resolution is computed in two steps:                                %    
% 1. Convolution in the modal space                                       %
% 2. Translation of the modal solutions into the original space           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = zeros(14,length(t));
eta = zeros(14,length(t)*2-1);
for i = 1:14
    h(i,:) = exp(-D(i,i)*t);
    eta(i,:) = conv(Ft(i,:),h(i,:))*t(2);
end
eta = eta(:,1:length(t)); % this is the only part of eta which is phisically relevant
y = Vn*eta; 

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
plot(t,y(11,:),'b','Linewidth',1.5)
plot(t,y(12,:),'r','Linewidth',1.5)
plot(t,y(13,:),'g','Linewidth',1.5)
plot(t,y(14,:),'m','Linewidth',1.5)

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

% representation of the front wheels' displacement

plot(t,y(11,:),'b','Linewidth',1.5)
plot(t,y(12,:),'r','Linewidth',1.5)

% representation is in [mm] scale 
% the four wheels' displacement is represented

xlabel('Time [s]')
legend('y1','y2')
ylabel('Displacement [mm]')
title('Wheel displacement, front wheels')
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

% representation of the rear wheels' displacement

plot(t,y(13,:),'m','Linewidth',1.5)
plot(t,y(14,:),'g','Linewidth',1.5)

% representation is in [mm] scale 
% the four wheels' displacement is represented

xlabel('Time [s]')
legend('y3','y4')
ylabel('Displacement [mm]')
title('Wheel displacement, rear wheels')
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

plot(t,y(9,:)/1000*57.2958,'g','Linewidth',2)
plot(t,y(10,:)/1000*57.2958,'r','Linewidth',2)

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
dphi = diff(y(9,:))./diff(t);
dphi2 = diff(dphi)./diff(t(1:end-1));
dtheta = diff(y(10,:))./diff(t);
dtheta2 = diff(dtheta)./diff(t(1:end-1));

hold on
plot(t(1:end-2),dphi2/1000,'g','Linewidth',1.5)
plot(t(1:end-2),dtheta2/1000,'r','Linewidth',1.5)
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
plot(t,y(8,:),'m','Linewidth',1.5)

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

subplot(2,1,2);
dz = diff(y(8,:))./diff(t);
dz2 = diff(dz)./diff(t(1:end-1));

hold on
plot(t(1:end-2),dz2/1000,'m','Linewidth',1.5)
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