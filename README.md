# NUM-Method-Project
clear, clc, close all
%%%%%%%%%% Driver Script
% Initial Parameters
load EnvironmentalForcing.mat
h = 0.01;
beta = 1;
mu_L = 6; 
mu_I = 10; 
e = 0.001;
A_p = 5000;
P_i = 930.27249;
S_i = P_i/A_p;
L_i = 0.01*S_i;
I_i = 0;
R_i = mu_I*I_i;
B_i = 1;
y0 = [S_i,L_i,I_i,R_i,B_i,0];
T_E = -0.35968 + 0.10789*T - 0.00214*(T.^2);
% Time span
tspan = 0:60;
[t, y] = rk4(@odefun, tspan, y0);



%%%%%%%%%% Functions

% ODEFUN function
function dydt = odefun(t, y)
S = y(1);
I = y(2);
L = y(3);
R = y(4);
P = y(5);
P_b = y(6);

% System of ODEs
dSdt = -beta*S*I + dPdt*(1/A_p);
dLdt = beta*S*I-(1/mu_L)*L + e;
dIdt = (1/mu_L) - (1/mu_I)*I;
dRdt = (1/mu_I)*I;
dP_bdt = (0.1724*P_b - 0.0000212*(P_b^2))*T_E;
dPdt = dP_bdt + (1.33*t_d)*T_E;
dydt = [dSdt; dIdt; dLdt; dRdt; dPdt; dP_bdt];
end

% SLIRP FUNCTION
function slopes = SLIRP

end

% RUNGE KUTTA 4TH ORDER FUNCTION
function [t, y] = rk4(odefun, tspan, y0)
 t0 = tspan(1);
tf = tspan(end);
h = tspan(2) - tspan(1);
num_steps = length(tspan);

t = zeros(num_steps, 1);
y = zeros(num_steps, length(y0));

% initial Conditions
t(1) = t0;
y(1, :) = y0;

% Runge Kutta
    for i = 2:num_steps
        t(i) = t(i-1) + h;
        k1 = h*odefun(t(i-1), y(i-1, :));
        k2 = h*odefun(t(i-1) + 0.5*h, y(i-1, :) + 0.5*k1);
        k3 = h*odefun(t(i-1) + 0.5*h, y(i-1, :)+0.5*k2);
        k4 = h*odefun(t(i-1) + h, y(i-1, :) + k3);
        
        y(i, :) = y(i-1, :)+(1/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%NEW AND IMPROVED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Driver Script
load EnvironmentalForcing.mat
beta_max = 1;
mu_L_min = 6; 
mu_I = 10; 
e = 0.001;
A_p = 5000;
P_i = 930.27249;
S_i = P_i/A_p;
L_i = 0.01*S_i;
I_i = 0;
R_i = mu_I*I_i;
B_i = 1;
y = [S_i,L_i,I_i,R_i,B_i,0];
y0 = y;
t0 = 0;
t_d = 1:length(T);

T_B = zeros(size(T)); 
for i = 1:length(T)
    if T(i) <= 0
        T_B(i) = 0;
    elseif T(i) > 0 && T(i) < 35
        T_B(i) =  (0.00241.*(T(i).^2.06737)).*(35 - T(i)).^(0.72859);
    elseif T(i) >= 35
        T_B(i) = 0;
    end
end

beta = beta_max*T_B;

mu_L = 0;
while mu_L < mu_L_min
    for i = 1:length(t_d)
        tinit = 0;
mu_L(i) = sum(tinit:T_B(i));
    end
end


slopes = @(t, y) SLIRP(t, y, beta, mu_L, mu_I, e, T, t_d, A_p);
[t, y] = rk4(slopes, tspan, y0);


% Plotting

plot(t, y(1,:),'k-')
hold on
plot(t, y(2, :), 'b-')
hold on
plot(t, y(3, :), 'g-')
hold on
plot(t, y(4, :), 'r-')
hold on
plot(t, y(5, :), 'p-')
hold on
plot(t, y(6, :), 'm-')
hold off

%% Functions


function slopes = SLIRP(t, y, beta, mu_L, mu_I, e, T, t_d, A_p)

S = y(1);
I = y(3);
L = y(2);
R = y(4);
P = y(5);
P_b = y(6);

T_E = -0.35968 + 0.10789.*T - 0.00214.*(T.^2);
dP_ldt = (1.33.*t_d).*T_E;
dLdt = beta.*S.*I-(1./mu_L).*L + e;
dIdt = (1./mu_L) - (1./mu_I).*I;
dRdt = (1./mu_I).*I;
dP_bdt = (0.1724.*P_b - 0.0000212.*(P_b.^2)).*T_E;
dPdt = dP_bdt + dP_ldt;
dSdt = -beta.*S.*I + dPdt.*(1/A_p);

slopes = [dSdt, dLdt, dIdt, dRdt, dP_bdt, dPdt];
end


function [t,y] = rk4(odefun,tspan,y0)
%UNTITLED Summary of this function goes here
% Detailed explanation goes here
num_steps = length(tspan);
t = zeros(1, num_steps);
y = zeros(length(y0), num_steps);
y(:, 1) = y0;

% Loop through time steps
for i = 1:num_steps-1
    h = tspan(i+1) - tspan(i);
    k1 = odefun(tspan(i), y(:, i));
    k2 = odefun(tspan(i) + h/2, y(:, i) + h*k1/2);
    k3 = odefun(tspan(i) + h/2, y(:, i) + h*k2/2);
    k4 = odefun(tspan(i) + h, y(:, i) + k3*h);

    y(:, i+1) = y(:, i) + (k1(i) + 2*k2(i) + 2*k3(i) + k4(i))/6;
    t(i+1) = tspan(i+1);
end
end
