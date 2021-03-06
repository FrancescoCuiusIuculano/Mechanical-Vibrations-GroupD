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

U = 20.3*ones(1,points+1);  % Nominal fluid velocity [m/s]

%% Definition of the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A pinned-pinned beam is considered: Galerkin method is used, with sine  %
% series as trial function. [w(x,t) = sum(q_i(t)*sin(i*pi*x))]            %
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
phi = zeros(N,length(x));
for i = 1:N
    phi(i,:) = sin(i*pi*x/L);
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

  % Point force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An impact force applied in the mid-point of the pipe is considered in   %
% this test. The force is modeled as a short sine half-wave.              %
% Such smoothness is required to avoid numerical errors with ode45.       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms t
F = sym(zeros(N,1));
syms sinbump(W,H,P)
sinbump(W,H,P) = H*(0.5+0.5*cos(2*pi*((t-P)/W)))*((P-W/2 < t) & (t < P+W/2));
f = sinbump(0.01,10,1);
for i = 1:N
    F(i) = -f*phi(i,fix(length(phi)/2));
end
Fp = [F; zeros(N,1)];

%% Resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The discrete equations are defined symbolically, and then they are      %
% converted into a vector field and resolved with Matlab's ode45 solver.  %
% The equilibrium position is used as initial condition.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms q(t) [N 1]                     % Definition of the unknowns
p(t) = [diff(q); q];                % Definition of the unknowns in state-space

eq = diff(p) == A*p + Fp;           % Definition of the system of equations

[Z,S] = odeToVectorField(eq);       % Conversion into vector field
% Z is the actual equation to solve, while S is a vector containing the
% substitutions made by odeToVectorField (which are not trivial). It will
% be useful later, to know which components of the solution to plot.
FUN = matlabFunction(Z,'vars', {'t','Y'}); % Conversion into function handle

interval = linspace(0,10,points);   % Time interval considered
yInit = zeros(1,2*N);               % Initial conditions
opt = odeset('Refine',15,'RelTol',1e-10,'AbsTol',1e-10,'InitialStep',1e-3,'MaxStep',1e-3); % To reduce numerical drift
pSol = ode45(FUN,interval,yInit,opt);

q1_val = deval(pSol,interval,3);
q2_val = deval(pSol,interval,1);
q3_val = deval(pSol,interval,5);
q4_val = deval(pSol,interval,7);
q5_val = deval(pSol,interval,9);
q6_val = deval(pSol,interval,11);
q7_val = deval(pSol,interval,13);

% Animating the response of the system
%%% WARNING: since the plot is updated repeatedly inside the for loop,
%%% closing the picture before the end of the run will make it re-open
%%% immediately. To stop the script, use ^C
fig = figure;
ff = matlabFunction(f);             % Conversion of the force into matlabFunction, to show the arrow
z = ones(1,points);

for tt = 1:5:length(interval)
    w = q1_val(tt)*phi(1,:) + q2_val(tt)*phi(2,:) + q3_val(tt)*phi(3,:) + ...
        q4_val(tt)*phi(4,:) + q5_val(tt)*phi(5,:) + q6_val(tt)*phi(6,:) + q7_val(tt)*phi(7,:);
    %plot3(z,x(1:length(w)).*1000,w*1e3,'Linewidth',1.5);
    tubeplot([z;x(1:length(w)).*1000;w*1e3*10],9.5,20);
%     daspect([0.6,5,0.4]);
%     colormap winter;
%     if interval(tt) < 1
%         hold on
%         quiver3(2+(1-interval(tt)),L/2,z,0,-2,0,'Linewidth',1.5) % Show a vertical arrow in the point where the force is applied
%         hold off
%     end
    
    xlim([-120 120])
    zlim([-120 120])
    ylim([0 2000])
    grid on
    title('Response to a pulse force in the mid-point')
    xlabel('width [mm]')
    ylabel('length [mm]')
    zlabel('height [m*10^-4]')
    drawnow
%     frame = getframe(fig);
%     im{tt}=frame2im(frame);
%     [AA,BB]=rgb2ind(im{tt},256);
%     if tt==1
%         imwrite(AA,BB,'ciao.gif','gif','LoopCount',Inf);
%     else
%         imwrite(AA,BB,'ciao.gif','gif','WriteMode','append');
%     end
end

%Tube representation -> function used
%https://it.mathworks.com/matlabcentral/fileexchange/5562-tubeplot

function [x,y,z]=tubeplot(curve,r,n,ct)
% Usage: [x,y,z]=tubeplot(curve,r,n,ct)
% 
% Tubeplot constructs a tube, or warped cylinder, along
% any 3D curve, much like the build in cylinder function.
% If no output are requested, the tube is plotted.
% Otherwise, you can plot by using surf(x,y,z);
%
% Example of use:
% t=linspace(0,2*pi,50);
% tubeplot([cos(t);sin(t);0.2*(t-pi).^2],0.1);
% daspect([1,1,1]); camlight;
%
% Arguments:
% curve: [3,N] vector of curve data
% r      the radius of the tube
% n      number of points to use on circumference. Defaults to 8
% ct     threshold for collapsing points. Defaults to r/2 
%
% The algorithms fails if you have bends beyond 90 degrees.
% Janus H. Wesenberg, july 2004
  if nargin<3 || isempty(n), n=8;
     if nargin<2, error('Give at least curve and radius');
     end
  end
  if size(curve,1)~=3
    error('Malformed curve: should be [3,N]');
  end
  if nargin<4 || isempty(ct)
    ct=0.5*r;
  end
  
  %Collapse points within 0.5 r of each other
  npoints=1;
  for k=2:(size(curve,2)-1)
    if norm(curve(:,k)-curve(:,npoints))>ct
      npoints=npoints+1;
      curve(:,npoints)=curve(:,k);
    end
  end
  %Always include endpoint
  if norm(curve(:,end)-curve(:,npoints))>0
    npoints=npoints+1;
    curve(:,npoints)=curve(:,end);
  end
  %deltavecs: average for internal points.
  %           first strecth for endpoitns.
  dv=curve(:,[2:end,end])-curve(:,[1,1:end-1]);
  %make nvec not parallel to dv(:,1)
  nvec=zeros(3,1);
  [buf,idx]=min(abs(dv(:,1))); nvec(idx)=1;
  xyz=repmat([0],[3,n+1,npoints+2]);
  
  %precalculate cos and sing factors:
  cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
  sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);
  
  %Main loop: propagate the normal (nvec) along the tube
  for k=1:npoints
    convec=cross(nvec,dv(:,k));
    convec=convec./norm(convec);
    nvec=cross(dv(:,k),convec);
    nvec=nvec./norm(nvec);
    %update xyz:
    xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1])+...
        cfact.*repmat(r*nvec,[1,n+1])...
        +sfact.*repmat(r*convec,[1,n+1]);
  end
  
  %finally, cap the ends:
  xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
  xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);
  
  %,extract results:
  x=squeeze(xyz(1,:,:));
  y=squeeze(xyz(2,:,:));
  z=squeeze(xyz(3,:,:));
  
  %... and plot:
  
  if nargout<3
      surf(x,y,z); 
      set(surf(x,y,z),'LineStyle','none')
  end
end