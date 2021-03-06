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
points = 10000;             % Number of discretization points

%% Definition of the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A pinned-pinned beam is considered: Galerkin method is used, with sine  %
% series as trial function. [w(x,t) = sum(q_i(t)*sin(i*pi*x))]            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = linspace(0,L,points);   % Longitudinal coordinate on the pipe
% Since Matlab's diff function used for the calculation of derivatives
% reduces the length of the vector by one, some dummy terms are added at
% the end of the vector x. These terms have no phisical meaning and will be
% removed later. Since fourth derivative is the highest derivative needed,
% four terms are added
for i = 1:4
    x(end+1) = x(end)+x(2);
end

%% Definition of the discrete problem
% Trial function
phi = zeros(N,length(x));
for i = 1:N
    phi(i,:) = sin(i*pi*x/L);
end
Dphi = diff(phi,1,2)/x(2);          % First derivative
Dphi2 = diff(phi,2,2)/x(2)^2;       % Second derivative
Dphi4 = diff(phi,4,2)/x(2)^4;       % Fourth derivative
% Vectors exceeding the length of interest are now rescaled
phi = phi(:,1:points);
Dphi = Dphi(:,1:points);
Dphi2 = Dphi2(:,1:points);

% Mass matrix
M = zeros(N);
for i = 1:N
    for j = 1:N
        M(i,j) = trapz(x(2),(mf+m)*phi(j,:).*phi(i,:));
    end
end

%% Resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In order to evaluate the behavior of the system with changing fluid     %
% velocity, a for cycle is used and the parameters depending by it are    %
% re-calculated at every iteration. This may be slow if high accuracy is  %
% required.                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num = 0;                            % Counter to keep track of the iterations
omega_n = zeros(N,5001);            % Natural circular frequency
for U = 0:.01:50
    num = num + 1;                  % Counter update
    % Calculation of the Reynolds number
    Re = U*Di/ni;
    if Re == 0
        f = 0;
    elseif Re < 2300
        f = 75/Re;
    else
        f = 0.3164/(Re^0.25);
    end
    % Stiffness Matrix
    K = zeros(N);
    for i = 1:N
        for j = 1:N
              K(i,j) = trapz(x(2),(E*I*Dphi4(j,:) + (mf*U^2*(f*L/(4*Di) + 1) - T + P*A*(1-2*v) )*Dphi2(j,:)).*phi(i,:));
        end
    end

    lambda = eig(M\K);              % Calculation of eigenvalues
    lambda = sort(lambda);
    omega_n(:,num) = sqrt(lambda);  % Calculation of natural circular frequencies
end

freq = omega_n/(2*pi);              % [rad/s] --> [Hz]

%% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The lowest natural frequency is represented as a function of the fluid's%
% velocity. A vertical line is used to highlight the critical velocity.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = linspace(0,U,num);              % Definition of the range of velocity
plot(u,freq(1,:),'Linewidth',1.5);  % Plot of f_1 vs U

% Search for the critical velocity
for i = 1:length(u)
    if real(freq(1,i)) == 0 && real(freq(1,i-1)) > 0
        U_crit = u(i);              % Save the value of the critical velocity 
        break                       % Quit the cycle
    end
end
hold on
xline(U_crit,'--k')                 % Draw the line for critical velocity
title('f_1 vs. U')
xlabel('U [m/s]')
ylabel('f_1 [Hz]')
text(U_crit,freq(1,1),strcat(' U_c =  ',num2str(U_crit),' m/s'),'VerticalAlignment','top','HorizontalAlignment','left')
