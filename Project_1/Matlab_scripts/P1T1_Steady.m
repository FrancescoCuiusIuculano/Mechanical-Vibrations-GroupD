close all;clear;clc;
%% Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The engine considered in this example is the Honda f20c (ap1),          %
% a 4 cylinders  inline engine with a 1,997cc displacemet, which equipped %
% the honda S2000 from 1999 to 2004.                                      %
% Source: https://en.wikipedia.org/wiki/Honda_F20C_engine                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Masses, source: https://www.s2ki.com/forums/s2000-under-hood-22/stock-f20c-f22c-piston-rod-weights-776633/
M = 148;                % Mass of the engine (kg), source: https://www.s2ki.com/forums/s2000-talk-1/how-much-does-f20c-weigh-64116/
mp = 0.355;             % Piston (kg)
mw = 0.109;             % Wrist pin (kg)
mr = 0.636;             % Rod (kg)
m = mp+mw+(2/3)*mr;     % Total oscillating mass
% Geometry
s = 0.084;              % Stroke (m)
r = s/2;
l = 0.153;              % Rod length (m)

%% System characteristics
k = 150000;             % Stiffness (N/m)
c = 800;                % Damping (Ns/m)
ni = 800;               % rpm in idle
nf = 8900;              % Limiter rpm

omega_n = sqrt(k/M);    % Natural frequency
zeta = c/(2*sqrt(k*M)); % Damping ratio
omega_d = omega_n*sqrt(1-zeta^2);

%% Sampling
Fs = 2048;              % Sampling frequency
dt = 1/Fs;              % Sampling period
L = 5;                  % Length of signal
t = 0:dt:L;             % Time vector

%% Force
omega_i = ni*2*pi/60;   % Rotational speed in idle (rad/s)
omega_f = nf*2*pi/60;   % Rotational speed at redline (rad/s)
rev = omega_i;          % In this test the engine is idling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In a 4 cylinders engine the first order inertial forces are always      %
% balanced, the only unbalanced forces are the second order ones          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definition of the forcing in terms of symbolic variables: this is needed 
% to solve our system's motion equation using Matlab's dsolve function
syms t_sym                                  % Time symbolic variable 
fi_sym = m*rev^2*r^2/l*cos(2*rev*t_sym);    % Second order inertial force
f_sym = 4*fi_sym;                           % Because there are 4 cylinders

% Definition of the forcing in terms of discrete variables: this is needed
% to compute the convolution using Matlab's conv function
fi = m*rev^2*r^2/l*cos(2*rev*t);            % Second order inertial force
f = 4*fi;                                   % Because there are 4 cylinders

%% Response analysis #1: Symbolic ode resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symbolic resolution of the equation of motion is computationally        %
% expensive. However, it allows to find a reliable reference in order to  %
% validate the results of the numerical approach.                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms y_sym(t_sym)                                                       % Displacement symbolic variable
Dy=diff(y_sym);                                                         % Derivative of y_sym
cond1 = y_sym(0)==0; cond2 = Dy(0)==0;                                  % Boundary conditions
conds=[cond1 cond2];                                            
ode = M*diff(y_sym,t_sym,2) + c*diff(y_sym,t_sym) + k*y_sym == f_sym;   % Equation of motion
ySol(t_sym) = dsolve(ode,conds);                                        % Resolution of the ode
ySol = simplify(ySol);                                                  % To speed up computations for plotting

%% Response analysis #2: Convolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical resolution using conv is much faster to compute than dsolve,  %
% but it requires the knowledge of the function h.                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = 1/(M*omega_d)*exp(-zeta*omega_n*t).*sin(omega_d*t); % Response to unitary impulse
y = conv(f,h)*t(2);                                     % Convolution ==> Response of the system

%% Plot

fs = 12;      % Fontsize
lw = 1;       % LineWidth
ms = 10;      % MarkerSize

%% Force plot representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In Figure 1, the force acting on the system is represented. The graph   %
% at the top shows the data in the entire interval taken into account,    %
% while the graph at the bottom contains a zoomed-in view of the same     %
% information, for better clarity.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(2,1,1)

plot(t,f/1000, 'LineWidth',0.1);
grid on
box on
title('2nd order excitation Force');
xlabel('Time [s]');
ylabel('Force [kN]');
set(gca, 'FontSize', fs, 'LineWidth', lw);

%% Force plot representation rescaled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2)

plot(t(1:fix(end/5)),f(1:fix(end/5))/1000,'-m', 'LineWidth',0.1);
grid on
box on
title('2nd order excitation Force (rescaled in time)');
xlabel('Time [s]');
ylabel('Force [kN]');
set(gca, 'FontSize', fs, 'LineWidth', lw);
set(gcf, 'Position', get(0, 'Screensize'));

%% Displacement using the symbolic approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In Figure 2, the results of the computations for the system response are%
% shown, toghether with a comparison of the results of both methods.      %
% As for figure 1, zoomed-in views are provided for better clarity.       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
subplot(3,2,1)

plot(t,ySol(t)*1000)
grid on
box on
title({'Displacement in time-domain obtainted using','the symbolic approach'});
xlabel('Time [s]');
ylabel('Displacement [mm]');
set(gca, 'FontSize', fs, 'LineWidth', lw);
set(gcf, 'Position', get(0, 'Screensize'));
%% Displacement using the symbolic approach rescaled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,2)

plot(t(1:fix(end/5)),ySol(t(1:fix(length(t)/5)))*1000,'-m');
grid on
box on
title({'Displacement in time-domain obtainted using','the symbolic approach (rescaled)'});
xlabel('Time [s]');
ylabel('Displacement [mm]');
set(gca, 'FontSize', fs, 'LineWidth', lw);
set(gcf, 'Position', get(0, 'Screensize'));

%% Displacement using the convolution approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,3)

plot(t,y(1:length(t))*1000)
grid on
box on
title({'Displacement in time-domain obtainted using','the convolution approach'});
xlabel('Time [s]');
ylabel('Displacement [mm]');
set(gca, 'FontSize', fs, 'LineWidth', lw);
set(gcf, 'Position', get(0, 'Screensize'));
%% Displacement using the convolution approach rescaled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,4)

plot(t(1:fix(end/5)),y(1:fix(length(t)/5))*1000,'-m');
grid on
box on
title({'Displacement in time-domain obtainted using','the convolution approach (rescaled)'});
xlabel('Time [s]');
ylabel('Displacement [mm]');
set(gca, 'FontSize', fs, 'LineWidth', lw);
set(gcf, 'Position', get(0, 'Screensize'));

%% Comparison between the two methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,5)

plot(t,ySol(t)*1000)
hold on
plot(t,y(1:length(t))*1000,'--','MarkerIndices',1:ms*2:length(t))
grid on
box on
title({'Displacement in time-domain:','dsolve vs conv'});
xlabel('Time [s]');
ylabel('Displacement [mm]');
legend('dsolve','conv')
set(gca, 'FontSize', fs, 'LineWidth', lw);
set(gcf, 'Position', get(0, 'Screensize'));
%% Comparison between the two methods rescaled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,6)

plot(t(fix(end-end/15):fix(end)),y(fix(length(t)-length(t)/15):fix(length(t)))*1000,'-m')%,'LineWidth',1.5);
hold on
plot(t(fix(end-end/15):fix(end)),ySol(t(fix(end-end/15):fix(end)))*1000,'og','MarkerIndices',1:ms:length(t));
grid on
box on
title({'Displacement in time-domain:','dsolve vs conv (rescaled, steady state)'});
xlabel('Time [s]');
ylabel('Displacement [mm]');
legend('dsolve','conv')
set(gca, 'FontSize', fs, 'LineWidth', lw);
set(gcf, 'Position', get(0, 'Screensize'));


%% Transmissibility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3 shows a plot of the transmissibility as a function of the      %
% circular frequency, highlighting the point in which the system is when  %
% the engine is idling.                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
tau = transmissibility(omega_n,zeta,omega_f,dt*100);
hold on
plot((2*rev/omega_n)^2,sqrt((2*zeta*2*rev/omega_n)^2+1)/sqrt((1-(2*rev/omega_n)^2)^2+(2*zeta*2*rev/omega_n)^2),'*r','LineWidth',2)
text((2*rev/omega_n)^2,sqrt((2*zeta*2*rev/omega_n)^2+1)/sqrt((1-(2*rev/omega_n)^2)^2+(2*zeta*2*rev/omega_n)^2),'(2\Omega/\omega_n)^2','VerticalAlignment','bottom','HorizontalAlignment','left')
grid on
box on
title({'Transmissibility'});
xlabel('(\omega/\omega_n)^2');
xlim([0 50])
ylabel('\tau (transmissibility)');
set(gca, 'FontSize', fs, 'LineWidth', lw);
set(gcf, 'Position', get(0, 'Screensize'));


function tau=transmissibility(omega_n, zeta, omega_f, resolution)
    omega = 0:resolution:fix(omega_f);
    tau = sqrt(((2*zeta*omega)./omega_n).^2+1)./sqrt((1-(omega./omega_n).^2).^2+(2*zeta*(omega./omega_n)).^2);
    plot((omega./omega_n).^2,tau,'LineWidth',1.5);
    hold on
    plot((omega./omega_n).^2,ones(1,length(omega)),':m','LineWidth',1.5)
end