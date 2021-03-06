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

N  = 7;                     % Number of terms to use

%% Definition of the problem
syms x t U 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A pinned-pinned beam is considered: Galerkin method is used, with sine  %
% series as trial function. [w(x,t) = sum(q_i(t)*sin(i*pi*x))]            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trial function
phi = sym(zeros(N,1));
for i = 1:N
    phi(i) = sin(i*pi*x/L);
end
% Mass matrix
M = sym(zeros(N));
for i = 1:N
    for j = 1:N
        M(i,j) = int((mf+m)*phi(j)*phi(i),x,0,L);
    end
end
% Stiffness matrix
% Calculation of the Reynolds number
Re = U*Di/ni;
f = 0.3164/(Re^0.25);
K = sym(zeros(N));
for i = 1:N
    for j = 1:N
        K(i,j) = int((E*I*diff(phi(j),x,4) + (mf*U^2*(f*L/(4*Di) + 1) - T + P*A*(1-2*v) )*diff(phi(j),x,2))*phi(i),x,0,L);
    end
end
% Damping matrix
C = sym(zeros(N));
for i = 1:N
    for j = 1:N
        C(i,j) = int(2*mf*U*diff(phi(j),x,1)*phi(i),x,0,L);
    end
end        

v = 20.3;                           % Nominal fluid velocity [m/s]
K = subs(K,U,v);  
C = subs(C,U,v);

% State-space matrix
A = [-M\C    -M\K;
      eye(N)  zeros(N)];
    
AA = double(A);                      % Conversion to double to improve performance
[V, D] = eig(AA);                    % Calculation of eigenvectors and eigenvalues
% Sorting the eigenvectors and the eigenvalues
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);
for i = 1:2*N
    eig = 1/Vs(i,i);
    for j = 1:2*N
        Vs(j,i) = Vs(j,i)*eig;
    end
end

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

[Z,S] = odeToVectorField(eq);       % Conversion into vector field
% Z is the actual equation to solve, while S is a vector containing the
% substitutions made by odeToVectorField (which are not trivial). It will
% be useful later, to know which components of the solution to plot.
FUN = matlabFunction(Z,'vars', {'t','Y'}); % Conversion into function handle

interval = linspace(0,10,10000);    % Time interval considered
mode = 1;                           % Desired mode (odd numbers: 1 for first mode, 3 for second ...)
yInit = [Vs(9,mode); Vs(2,mode);  Vs(8,mode); Vs(1,mode);  Vs(11,mode); Vs(4,mode);  Vs(13,mode);
         Vs(6,mode); Vs(10,mode); Vs(3,mode); Vs(12,mode); Vs(5,mode);  Vs(14,mode); Vs(7,mode)]; % Initial conditions
opt = odeset('Refine',15,'RelTol',1e-10,'AbsTol',1e-10,'InitialStep',1e-3,'MaxStep',1e-2); % To reduce numerical drift
pSol = ode45(FUN,interval,yInit,opt);

q1_val = deval(pSol,interval,3);
q2_val = deval(pSol,interval,1);
q3_val = deval(pSol,interval,9);
q4_val = deval(pSol,interval,5);
q5_val = deval(pSol,interval,11);
q6_val = deval(pSol,interval,7);
q7_val = deval(pSol,interval,13);

% Plotting of the desired mode
pipe = linspace(0,L,10000);         % Vector representing the pipe
phi = double(subs(phi,x,pipe));     % The conversion of phi to double is done outside of the loop to improve performance

hold on

colors = ["#E3F2FD"; "#BBDEFB"; "#90CAF9"; "#64B5F6"; "#42A5F5";
          "#2196F3"; "#1E88E5"; "#1976D2"; "#1565C0"; "#0D47A1"];
cc = 0;
for tt = 1:1000:length(interval)
    cc = cc + 1;
    w = q1_val(tt)*phi(1,:) + q2_val(tt)*phi(2,:) + q3_val(tt)*phi(3,:) + ...
        q4_val(tt)*phi(4,:) + q5_val(tt)*phi(5,:) + q6_val(tt)*phi(6,:) + q7_val(tt)*phi(7,:);
    plot(pipe,w*1e3,'color',colors(cc),'Linewidth',1.5);
end
title('1st mode in t (9sec)')
xlabel('x [m]');
ylabel('w [mm]');
legend('t = 0 s','t = 1 s','t = 2 s','t = 3 s','t = 4 s','t = 5 s','t = 6 s','t = 7 s','t = 8 s','t = 9 s')