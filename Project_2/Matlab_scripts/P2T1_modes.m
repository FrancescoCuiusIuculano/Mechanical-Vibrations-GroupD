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

%% Analysis of the undamped system

%When dealing with undamped systems the natural frequencies can be found
%straigth forward, by determining the eigen values of the system. The modes
%of vibration can be simlutaneously calculated as part of the solving
%method.

[U,lambda] = eig(M\K);                                    % Calculation of eigenvectors and eigenvalues
 Mn = transpose(U)*M*U;                                % Matrix that contains in its diagonal the modal masses
 
% Normalizing the modal matrix with the help of the modal mass
Un = zeros(7);                                                   
for j =1:7
    nom = 1/sqrt(Mn(j,j));
    for i = 1:7
        Un(i,j)=U(i,j)*nom;
    end
end

% Calculation of natural frequencies
omega = zeros(1,7);
for i = 1:7
    omega(i) = sqrt(lambda(i,i));                       %Creating a vector that contains the natural frequencies of the undamped 
end                                                                     %system

% Calculation of mode shapes
syms t                                                                %Creating a symbolic variable t for calculate the response as a function of t
q = sym(zeros(7));                                            % Pre-allocationg for speeding up the code
for i = 1:7
    q(:,i) = Un(:,i)*cos(omega(i)*t);                   % Generating the matrix that contains the response of the system 
end

time = 0:.001:1;
mode = zeros(7,length(time),7);                     % Pre-allocationg for speeding up the code
for j = 1:7
    for i = 1:7
        mode(i,:,j) = double(subs(q(i,j),t,time));  %Transform the symbolic variable into a vector and reating a matrix 
    end                                                                 %with respect to each mode of vibration that contains the response
end                                                                     % of every DOF with respect to each mode

% Plotting of mode shapes
for j = 1:7                                                           %For loop design to plot every DOF with respect to every mode
    subplot(4,2,j)
    hold on
    for i = 1:7
        plot(time,mode(i,:,j));
    end
    xlabel('Time (s)')
    ylabel('Displacement')
    title(strcat('Mode ', num2str(j)))
    grid minor
    set(gca, 'FontName', 'Helvetica')
    set(gca, 'FontSize',13)
    set(gca, 'GridAlpha',1)
    grid on
    set(gca, 'MinorGridAlpha',0.9)
    set(gcf, 'Position', get(0, 'Screensize')-[0 0 10 10]);
end
subplot(4,2,8)
plot(0,0,  0,0,  0,0,  0,0, 0,0, 0,0, 0,0)
axis off
legend({'z [m]','\phi [rad]','\theta [rad]','y1 [m]','y2 [m]','y3 [m]','y4 [m]'},'FontSize',13)

% Exporting the data for use in Blender
% Commented out to avoid spamming professors' PCs with files
% for j = 1:7
%     csvwrite(strcat(num2str(j), '_z.csv'),transpose(mode(1,:,j)))
%     csvwrite(strcat(num2str(j), '_y1.csv'),transpose(mode(4,:,j)))
%     csvwrite(strcat(num2str(j), '_y2.csv'),transpose(mode(5,:,j)))
%     csvwrite(strcat(num2str(j), '_y3.csv'),transpose(mode(6,:,j)))
%     csvwrite(strcat(num2str(j), '_y4.csv'),transpose(mode(7,:,j)))
%     csvwrite(strcat(num2str(j), '_angles.csv'),[transpose(mode(2,:,j)) transpose(mode(3,:,j))])
% end

%% Analysis of the damped system 
% This is for calculating the eigen vectors as well as the modes
% of vibration for the damped system
A = [zeros(7,7), M;
        M, C];
    
B = [-M, zeros(7,7);
        zeros(7,7), K];

[V, D] = eig(A\B);                                              % calculation of eigenvectors and eigenvalues
Ma= transpose(V)*A*V;                                   % matrix containing the complex modal masses in its diagonal 

% Normalizing the modal matrix
Vn= zeros([14, 14], 'like', Ma);                        % preallocating for optimizing speed
for j =1:14
    nom = 1/sqrt(Ma(j,j));
    for i = 1:14 
        Vn(i,j)=V(i,j)*nom;
    end
end

omega_d = zeros(1,14);                                      %Defining a vector that contains the complex eigenvalues
for i = 1:14
    omega_d(i) = abs(D(i,i));                                 %Filling up the vector with the complex eigenvalues
end

syms t                                                                %Setting the symbolic variable t in order to create a matrix y
y = sym(zeros(14,7));                                       % that contains the response of the system 
for i = 1:2:14
    y(:,ceil(i/2)) = 2*real(V(:,i)*exp(-D(i,i)*t));    % Defining the response with the help of the complex eigenvectors 
end                                                                     % as the equation (2.148) suggested (see the lecture notes for more
                                                                            %explanation)
qd = y(8:14,:);
shape = zeros(7,length(time),7);
for j = 1:7
    for i = 1:7
        shape(i,:,j) = double(subs(qd(i,j),t,time));%Changing the symbolic variable into a time vector and creating a matrix 
    end                                                                  % for obtaining the response of each DOF with respect to every single mode
end

figure(2)

for j = 1:7                                                           %For loop design to plot every DOF with respect to each complex eigenvector 
    subplot(4,2,j)
    hold on
    for i = 1:7
        plot(time,shape(i,:,j));
    end
    xlabel('Time (s)')
    ylabel('Displacement')
    title(strcat('Mode ', num2str(j)))
    grid minor
    set(gca, 'FontName', 'Helvetica')
    set(gca, 'FontSize',13)
    set(gca, 'GridAlpha',1)
    grid on
    set(gca, 'MinorGridAlpha',0.9)
    set(gcf, 'Position', get(0, 'Screensize')-[0 0 10 10]);
end
subplot(4,2,8)
plot(0,0,  0,0,  0,0,  0,0, 0,0, 0,0, 0,0)
axis off
legend({'z [m]','\phi [rad]','\theta [rad]','y1 [m]','y2 [m]','y3 [m]','y4 [m]'},'FontSize',13)


% Exporting the data for use in Blender
% Commented out to avoid spamming professors' PCs with files
% for j = 1:7
%     csvwrite(strcat(num2str(j), '_zd.csv'),transpose(shape(1,:,j)))
%     csvwrite(strcat(num2str(j), '_y1d.csv'),transpose(shape(4,:,j)))
%     csvwrite(strcat(num2str(j), '_y2d.csv'),transpose(shape(5,:,j)))
%     csvwrite(strcat(num2str(j), '_y3d.csv'),transpose(shape(6,:,j)))
%     csvwrite(strcat(num2str(j), '_y4d.csv'),transpose(shape(7,:,j)))
%     csvwrite(strcat(num2str(j), '_anglesd.csv'),[transpose(shape(2,:,j)) transpose(shape(3,:,j))])
% end
