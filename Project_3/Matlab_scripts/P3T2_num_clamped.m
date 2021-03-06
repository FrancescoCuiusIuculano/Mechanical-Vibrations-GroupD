clear;clc;close all;

%% Definition of data
% -1 assignment to enter the cicle. Commented value is the standard one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSERT THE MAIN PARAMETERS RESPECTING A GIVEN RANGE                     %
% If no value is inserted, the default one is used                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Pipe parameters input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Samco eXtreme SHL/19 is used as reference                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E  = -1;                    % Young modulus: 0.1e9 [Pa]
Do = -1;                    % Outer diameter: 19e-3; [m]
th = -1;                    % Wall thickness: 3.9e-3 [m]
L  = -1;                    % Length: 2; [m]
rho_p = -1;                 % Density: 2300; [kg/m^3]
v = -1;                     % Poisson's ratio: 0.49 [Dimensionless]

varSet=[0.1e9, 2, 19e-3, 3.9e-3, 2300, 0.49];
keys = {'E', 'L', 'Do', 'th', 'rho_p', 'v', 'E1', 'E2', 'L1', 'L2', 'Do1', 'Do2', 'th1', 'th2', 'rho_p1', 'rho_p2', 'v1', 'v2'};
values =[E L Do th rho_p  v 1e6 5e11 0.1 15 10e-3 1 1e-3 1e-2 1000 100000 0.01 0.5];
limits= containers.Map(keys, values, 'UniformValues', false);

fprintf('Leave blank to use the default value\n\n');
for e=1:6
  while ((limits(char(keys(e)))<=limits(char(strcat(keys(e), '1')))) || (limits(char(keys(e)))>limits(char(strcat(keys(e), '2')))))
     limits(char(keys(e))) = input(strcat("Insert the ", keys(e), " of the pipe (range ", string(limits(char(strcat(keys(e), '1')))), " - ", string(limits(char(strcat(keys(e), '2')))), ") = "));
    if isempty(limits(char(keys(e))))
        limits(char(keys(e)))=varSet(e);
    end
  end
end
E= limits('E');
Do = limits('Do');
th = limits('th');
L = limits('L');
rho_p = limits('rho_p');
v =limits('v');

Di = Do - 2*th;             % Inner diameter [m]
A = (Do^2-Di^2)*pi;         % Cross-section [m^2]
I = (pi/64)*(Do^4-Di^4);    % Moment of Inertia [m^4]
m = rho_p*pi*(Do^2-Di^2)/4; % Mass per unit of length [kg/m]

% Fluid parameters input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_f = -1;                 % Density: 1077; [kg/m^3]
ni = -1;                    % Kinematic viscosity: 3.98e-6 [m^2/s]

varSet_f = [1077, 3.98e-6];
keys_f = {'rho_f', 'ni',  'rho_f1', 'rho_f2', 'ni1','ni2'};
values_f = [rho_f ni 0.1 400000 1e-9 1e-1];
limits_f = containers.Map(keys_f, values_f, 'UniformValues', false);

for e=1:2
  while ((limits_f(char(keys_f(e)))<=limits_f(char(strcat(keys_f(e), '1')))) || (limits_f(char(keys_f(e)))>limits_f(char(strcat(keys_f(e), '2')))))
     limits_f(char(keys_f(e))) = input(strcat("Insert the ", keys_f(e), " of the fluid (range ", string(limits_f(char(strcat(keys_f(e), '1')))), " - ", string(limits_f(char(strcat(keys_f(e), '2')))), ") = "));
    if isempty(limits_f(char(keys_f(e))))
        limits_f(char(keys_f(e)))=varSet_f(e);
    end
  end
end
rho_f = limits_f('rho_f');
ni = limits_f('ni');

mf = rho_f*pi*Di^2/4;       % Mass per unit of length [kg/m]

% Other parameters
P = 2e5;                    % Pressure [Pa]
T = 100;                    % Tension [N]

N  = 4;                     % Number of terms to use
points = 10000;             % Number of discretization points

U = 20.3*ones(1,points+1);  % Nominal fluid velocity [m/s]

%% Definition of the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A clamped-clamped beam is considered: Galerkin method is used.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = linspace(0,L,points);       % Longitudinal coordinate on the pipe
% Since Matlab's diff function used for the calculation of derivatives
% reduces the length of the vector by one, some dummy terms are added at
% the end of the vector x. These terms have no phisical meaning and will be
% removed later. Since fourth derivative is the highest derivative needed,
% four terms are added
for i = 1:4
    x(end+1) = x(end)+x(2);
end

t = linspace(0,10,points);      % Time
t(end+1) = t(end) + t(2);       % Needeed to compensate diff, as above

U = U + 0.15*sin(100*2*pi*t);   % Velocity oscillation
DU = diff(U)/t(2);              % Derivative

% Trial function
beta = [4.730041 7.853205 10.995608 14.137165]/L;
alpha = (sinh(beta*L)-sin(beta*L))./(cos(beta*L)-cosh(beta*L));

phi = zeros(N,points+4);
for i = 1:N
    phi(i,:) = sinh(beta(i)*x) - sin(beta(i)*x) + alpha(i)*(cosh(beta(i)*x) - cos(beta(i)*x));
end
Dphi = diff(phi,1,2)/x(2);      % First derivative
Dphi2 = diff(phi,2,2)/x(2)^2;   % Second derivative
Dphi4 = diff(phi,4,2)/x(2)^4;   % Fourth derivative
% Vectors exceeding the length of interest are now rescaled
x = x(1:points);
phi = phi(:,1:points);
Dphi = Dphi(:,1:points);
Dphi2 = Dphi2(:,1:points);

U = U(1:points);

%% Definition of the discrete problem
% Mass matrix
M = zeros(N);
for i = 1:N
    for j = 1:N
        M(i,j) = trapz(x(2),(mf+m)*phi(j,:).*phi(i,:));
    end
end

% Calculation of the Reynolds number
Re = U*Di/ni;
if Re == 0
    f = 0;
elseif Re < 2300
    f = 75./Re;
else
    f = 0.3164./(Re.^0.25);
end
% Stiffness Matrix
K = zeros(N);
for i = 1:N
    for j = 1:N
        K(i,j) = trapz(x(2),(E*I*Dphi4(j,:) + (mf*U.^2.*(f*L/(4*Di) + 1) + DU.*(L-x) - T+P*A*(1-2*v)).*Dphi2(j,:)).*phi(i,:));
    end
end
% Damping matrix
C = zeros(N);
for i = 1:N
    for j = 1:N
        C(i,j) = trapz(x(2),2*mf*U.*Dphi(j,:).*phi(i,:));
    end
end

% State-space matrix
A = [-M\C    -M\K;
      eye(N)  zeros(N)];

[V, D] = eig(A);                    % Calculation of eigenvectors and eigenvalues
% Sorting the eigenvectors and the eigenvalues
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);

omega_d = zeros(1,2*N);             % Defining a vector that contains the complex eigenvalues
for i = 1:2*N
    omega_d(i) = abs(D(i,i));       % Filling up the vector with the complex eigenvalues
end

%% Resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The discrete equations are defined symbolically, and then they are      %
% converted into a vector field and resolved with Matlab's ode45 solver.  %
% To obtain the N-th mode, the N-th eigenvector is used as initial        %
% condition.                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms q(t) [N 1]                     % Definition of the unknowns
p(t) = [diff(q); q];                % Definition of the unknowns in state-space

eq = diff(p) == A*p;                % Definition of the system of equations

[Z,S] = odeToVectorField(eq); % Conversion into vector field
% Z is the actual equation to solve, while S is a vector containing the
% substitutions made by odeToVectorField (which are not trivial). It will
% be useful later, to know which components of the solution to plot.
FUN = matlabFunction(Z,'vars', {'t','Y'}); % Conversion into function handle

interval = linspace(0,15,points);   % Time interval considered
mode = 1;                           % Desired mode (odd numbers: 1 for first mode, 3 for second ...)
yInit = [Vs(6,mode); Vs(2,mode); Vs(5,mode); Vs(1,mode); Vs(7,mode); Vs(3,mode); Vs(8,mode); Vs(4,mode)];
opt = odeset('Refine',15,'RelTol',1e-10,'AbsTol',1e-10,'InitialStep',1e-3,'MaxStep',1e-2); % To reduce numerical drift
pSol = ode45(FUN,interval,yInit,opt);

q1_val = deval(pSol,interval,3);
q2_val = deval(pSol,interval,1);
q3_val = deval(pSol,interval,5);
q4_val = deval(pSol,interval,7);

% Plotting of the desired mode
hold on
colors = ["#E3F2FD"; "#BBDEFB"; "#90CAF9"; "#64B5F6"; "#42A5F5";
          "#2196F3"; "#1E88E5"; "#1976D2"; "#1565C0"; "#0D47A1"];
cc = 0;
for tt = 1:points/10:length(interval)
    cc = cc + 1;
    w = q1_val(tt)*phi(1,:) + q2_val(tt)*phi(2,:) + q3_val(tt)*phi(3,:) + q4_val(tt)*phi(4,:);
    plot(x(1:length(w)),w*1e3,'color',colors(cc),'Linewidth',1.5);
end
title('1st mode in t (13.5sec)')
xlabel('x [m]');
ylabel('w [mm]');
legend('t = 0 s','t = 1.5 s','t = 3 s','t = 4.5 s','t = 6 s','t = 7.5 s','t = 9 s','t = 10.5 s','t = 12 s','t = 13.5 s')