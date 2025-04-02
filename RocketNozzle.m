%% Coded for the fulfilment of Master's Degree at Politecnico Di Milano
% Author:: Eshwar Bhargav Bhupanam
% Course:: Modeling and Simulation of Aerospace Systems
% Topic:: Electrical Circuit - LCR
% Year:: 2020-2021

%% Initializing....
clear vars; clc; close all

%% Nozzle Thermal Test - Physical model inputs
par.A = 1;  % area of section m2
par.l = [0.005 0.03 0.025 0.005]; % section length in 'm'
par.k = [40 216 0.038 151]; % Thermal conductivity in 'W/mK' 1245
par.rho = [0 1850 105 0]; % Density of Section in 'kg/m3'   not to be considered becase capacitance can be negligible
par.c = [0 1825 835 0]; % Specific heat capacity of Section in 'J/kgK'

R = par.l./(par.k.*par.A);  % Section Resistance
C = par.rho.*par.A.*par.l.*par.c;   % Section Capacitance

T_cel = [20, 1000, 20, 20]; % [in_before firing, in_after firing, out_before firing, out_after firing]
[T_kel] = celsius_2_kelvin(T_cel);   % [in_before firing, in_after firing, out_before firing, out_after firing]

% Resistance between nodes
R_int = 0.01;   % Resistance for interface

% Case 1 (Single lumped case)
% Initial Conditions
tspan = linspace(0,60,1000);
% Variation of inner sheet temperature with time
[Ti] = temp_time_disc(tspan,T_kel,0);
% Variation of outer casing temperature with time
[To] = temp_time_disc(tspan,T_kel,1);

cond_deg = [20 20];
cond_kel = celsius_2_kelvin(cond_deg);
conds = cond_kel;
% Thermal system function
[t_single,T_single] = thermalsystem(R,C,R_int,T_kel,Ti,conds,tspan, 1);
delete Thermalodeset1.m
% Plotting - Temperature Profile - Single lumped case
figure()
hold on; grid on;
plot(t_single,[Ti T_single To]);
legend({'$T_{i}$','$T_{1}$','$T_{2}$','$T_{3}$','$T_{4}$','$T_{5}$','$T_{o}$'}, 'Interpreter', 'latex');
ylabel('Temperature [K]', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')

% Case 2 (Two lumped case)
% Initial conditions 
cond_deg = [20 20 20 20];
cond_kel = celsius_2_kelvin(cond_deg);
conds = cond_kel;
% Thermal system function
[t_two,T_two] = thermalsystem(R,C,R_int,T_kel,Ti,conds,tspan, 2);
delete Thermalodeset2.m
% Plotting - Temperature Profile - Two lumped case
figure(2)
plot(t_two,[Ti T_two To]);
legend({'$T_{i}$','$T_{1}$','$T_{2a}$','$T_{2b}$','$T_{3}$','$T_{4a}$','$T_{4b}$','$T_{5}$','$T_{o}$'}, 'Interpreter', 'latex');
ylabel('Temperature [K]', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')
grid on 

%% Functions

function [T_kelvin] = celsius_2_kelvin(T_celsius)
T_kelvin = 273+T_celsius;
end

function [T] = temp_time_disc(t,T_kelvin,type)
T = zeros(length(t),1);
for i=1:length(t)
    if type == 0 && t(i)<1
        T(i) = T_kelvin(1)+980*t(i);
    elseif type == 0 && t(i)>=1
        T(i,:) = T_kelvin(2);
    elseif type == 1
        T(i,:) = T_kelvin(4);
    end
end
end

function [t,T] = thermalsystem(R,C,R_int,T_kel,Ti,conds,tspan,type)
if type == 1
    R1 = R(1)+R(2)/2;
    R3 = R(2)/2+R_int+R(3)/2;
    R5 = R(3)/2+R(4);
    
    % Differential equations
    syms T2(t) T4(t)
    
    EOM1 = diff(T2) == 1/C(2)*((1/R1) *(T_kel(2)-T2) - (1/R3)*(T2-T4));
    EOM2 = diff(T4) == 1/C(3)*(((1/R3) *(T2-T4) - (1/R5)*(T4-T_kel(4))));
    
    EOM = [EOM1; EOM2];
    V = odeToVectorField(EOM);
    Thermalodeset1 = matlabFunction(V,'vars', {'t', 'Y'},'file','Thermalodeset1.m');
    
    % Integrating in the time domain
    options = odeset;
    [t,y] = ode45(@Thermalodeset1,tspan,conds,options);
    
    % Nodal temperatures
    T(:,2) = y(:,2);
    T(:,4) = y(:,1);
    
    %Heat Flow
    Qi2 = (1/(R(1)+R(2)/2))*(Ti - T(:,2));
    Q24 = (1/(R(2)/2+R_int+R(3)/2))*(T(:,2) - T(:,4));
    Q4o = (1/((R(4)/2))) *(T(:,4)-T_kel(4));
    
    % Nodal temperatures
    T(:,1) = Ti - Qi2*(R(1)/2);
    T(:,3) = T(:,2) - Q24*(R(2)/2+R_int/2);
    T(:,5) = T_kel(4) + Q4o*(R(4)/2);
elseif type == 2
    R1 = R(1)+R(2)/3;
    R2 = R(2)/3;
    R3 = R(2)/3+R_int+R(3)/3;
    R4 = R(3)/3;
    R5 = R(3)/3+R(4);
    
    syms T2a(t) T2b(t) T4a(t) T4b(t)
    
    EOM1 = diff(T2a) == 1/(C(2)/2)*((1/R1) *(T_kel(2)-T2a) - (1/R2)*(T2a-T2b));
    EOM2 = diff(T2b) == 1/(C(2)/2)*((1/R2) *(T2a-T2b) - (1/R3)*(T2b-T4a));
    EOM3 = diff(T4a) == 1/(C(3)/2)*((1/R3) *(T2b-T4a) - (1/R4)*(T4a-T4b));
    EOM4 = diff(T4b) ==1/(C(3)/2)*((1/R4) *(T4a-T4b) - (1/R5)*(T4b-T_kel(4)));
    
    EOM = [EOM1; EOM2; EOM3; EOM4];
    V = odeToVectorField(EOM);
    Thermalodeset2 = matlabFunction(V,'vars', {'t', 'Y'},'file','Thermalodeset2.m');
    
    % Integrating in the time domain
    options = odeset;
    [t,y] = ode45(@Thermalodeset2,tspan,conds,options);
    
    % Nodal Temperatures
    T(:,2) = y(:,2);
    T(:,3) = y(:,1);
    T(:,5) = y(:,3);
    T(:,6) = y(:,4);
    
    % Heat flow
    Qi2a = (1/(R(1)+R(2)/3))*(Ti - T(:,2));
    Q2b4a = (1/(R(2)/3+R_int+R(3)/3)) *(T(:,3) - T(:,5));
    Q4bo = (1/((R(4)/2))) *(T(:,6)-T_kel(4));
    
    % Nodal temperatures
    T(:,1) = Ti - Qi2a*(R(1)/2);
    T(:,4) = T(:,3) - Q2b4a*(R(2)/3+R_int/2);
    T(:,7) = T_kel(4) + Q4bo*(R(4)/3);
end
end
