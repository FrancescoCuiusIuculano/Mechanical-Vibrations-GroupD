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
rev = omega_f/L*t;      % In this test the engine is ramping form zero to redline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In a 4 cylinders engine the first order inertial forces are always      %
% balanced, the only unbalanced forces are the second order ones          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definition of the forcing in terms of discrete variables: this is needed
% to compute the convolution using Matlab's conv function
fi = m*rev.^2*r^2/l.*cos(2*rev.*t);         % Second order inertial force
f = 4*fi;                                   % Because there are 4 cylinders

%% Response analysis: Convolution
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

%% Displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In Figure 2, the results of the computations for the system response are%
% shown. As for figure 1, a zoomed-in view is provided for better clarity.%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
subplot(2,1,1)

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
subplot(2,1,2)

plot(t(1:fix(end/5)),y(1:fix(length(t)/5))*1000,'-m');
grid on
box on
title({'Displacement in time-domain obtainted using','the convolution approach (rescaled)'});
xlabel('Time [s]');
ylabel('Displacement [mm]');
set(gca, 'FontSize', fs, 'LineWidth', lw);
set(gcf, 'Position', get(0, 'Screensize'));

%% Transmissibility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3 shows a plot of the transmissibility as a function of the      %
% circular frequency, as well as the force transmitted to the chassis over%
% time.                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
subplot(2,1,1)
tau = transmissibility(omega_n,zeta,omega_f,omega_f/L*dt);
grid on
box on
title({'Transmissibility'});
xlabel('(\omega/\omega_{n})^2');
xlim([0 50])
ylabel('\tau (transmissibility)');
set(gca, 'FontSize', fs, 'LineWidth', lw);
set(gcf, 'Position', get(0, 'Screensize'));

subplot(2,1,2)
plot(t,4.*tau.*m.*rev.^2.*r.^2./l,'LineWidth',1.5);
grid on
box on
title({'Amplitude of transmitted force'});
xlabel('Time [s]');
ylabel('|T| [N]');
set(gca, 'FontSize', fs, 'LineWidth', lw);
set(gcf, 'Position', get(0, 'Screensize'));

function tau=transmissibility(omega_n, zeta, omega_f, resolution)
    omega = 0:resolution:omega_f;
    tau = sqrt(((2*zeta*omega)./omega_n).^2+1)./sqrt((1-(omega./omega_n).^2).^2+(2*zeta*(omega./omega_n)).^2);
    plot((omega./omega_n).^2,tau,'LineWidth',1.5);
    hold on
    plot((omega./omega_n).^2,ones(1,length(omega)),':m','LineWidth',1.5)
end