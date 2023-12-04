
% Sydney Nelson, Theodore Hays
% u1399500,      u1379892
% Lab 09
clear, clc, close all


load EnvironmentalForcing.mat

B_max = 1;
mu_L_min = 6; 
mu_I = 10; 
e = 0.001; 
A_P = 5000; 
P_I = 1.33*30*(-0.35968+0.10789*15-0.00214*15*15)*30; 
S_I = P_I/A_P; 
L_I = S_I*0.01;
I_I = 0; 
R_I = I_I*mu_I; 
B_I = 1; 

tinit=1;
T_B = zeros(size(T)); 
mu_L = zeros(size(T));

for i=1:length(T)

    T_B(i) = solve(T(i));
    mu_L(i) = sum(T_B(tinit:i));

    while(mu_L(i)>mu_L_min)

        tinit=tinit+1;
        mu_L(i) = sum(T_B(tinit:i));

    end
end
mu_L = 1./mu_L; 

p{1} = B_max; 
p{2} = mu_L; 
p{3} = mu_I; 
p{4} = e;    
p{5} = T;    
p{6} = tspan;
p{7} = A_P;    

odefun = @(t,y) slirp_function(t,y,p);

y0(1) = B_I; 
y0(2) = P_I; 
y0(3) = S_I; 
y0(4) = L_I;
y0(5) = I_I; 
y0(6) = R_I; 

[tspan,y] = rk4(odefun,tspan,y0');

B = y(:,1);
P = y(:,2);
S = y(:,3);
L = y(:,4);
I = y(:,5);
R = y(:,6);

hold on

plot(tspan,P/A_P,'-g',tspan,B/A_P,'--m',tspan,S,'-.k',tspan,L,'--c',tspan,I,':b',tspan,R,'-.r','LineWidth',2);
title('SLIRP model');
xlabel('days (t)');
ylabel('pop fraction')
legend('total pop','susceptible pop','berry pop','latent pop ','infected pop','removed pop','Location','northwest')
xlim([0,61]);
hold off

figure
Bvals = [1 1.5 2];
for k = 1:3
B_max = Bvals(k);
mu_L_min = 6; 
mu_I = 10; 
e = 0.001; 
A_P = 5000; 
P_I = 1.33*30*(-0.35968+0.10789*15-0.00214*15*15)*30; 
S_I = P_I/A_P; 
L_I = S_I*0.01;
I_I = 0; 
R_I = I_I*mu_I; 
B_I = 1; 

tinit=1;
T_B = zeros(size(T)); 
mu_L = zeros(size(T));
for i=1:length(T)
    T_B(i) = solve(T(i));
    mu_L(i) = sum(T_B(tinit:i));
    while(mu_L(i)>mu_L_min)
        tinit=tinit+1;
        mu_L(i) = sum(T_B(tinit:i));
    end
end
mu_L = 1./mu_L; 

p{1} = B_max; 
p{2} = mu_L; 
p{3} = mu_I; 
p{4} = e;    
p{5} = T;    
p{6} = tspan;
p{7} = A_P;

odefun = @(t,y) slirp_function(t,y,p);

y0(1) = B_I; 
y0(2) = P_I; 
y0(3) = S_I; 
y0(4) = L_I;
y0(5) = I_I; 
y0(6) = R_I; 
[tspan,y] = rk4(odefun,tspan,y0');
B = y(:,1);
P = y(:,2);
S = y(:,3);
L = y(:,4);
I = y(:,5);
R = y(:,6);

subplot(2, 3, k);
plot(tspan,P/A_P,'-g',tspan,B/A_P,'--m',tspan,S,'-.k',tspan,L,'--c',tspan,I,':b',tspan,R,'-.r','LineWidth',2);
title('Beta Max = ', B_max);
xlabel('days (t)');
ylabel('Population (fraction of initial)')
legend('total pop','susceptible pop','berry pop','latent pop ','infected pop','removed pop','Location','northwest')
xlim([0,61]);
hold on

end

mu_L_vals = [1 5 10];
for k = 1:3
P_I = 1.33*30*(-0.35968+0.10789*15-0.00214*15*15)*30; 
S_I = P_I/A_P; 
L_I = S_I*0.01;
I_I = 0; 
R_I = I_I*mu_I; 
B_I = 1; 
mu_I = 10; 
e = 0.001; 
mu_L_min = mu_L_vals(k);
B_max = 1;
A_P = 5000; 

tinit=1;
T_B= zeros(size(T)); 
mu_L= zeros(size(T));
for i=1:length(T)
    T_B(i) = solve(T(i));
    mu_L(i) = sum(T_B(tinit:i));
    while(mu_L(i)>mu_L_min)
        tinit=tinit+1;
        mu_L(i) = sum(T_B(tinit:i));
    end
end
mu_L = 1./mu_L; 

p{1} = B_max; 
p{2} = mu_L; 
p{3} = mu_I; 
p{4} = e;    
p{5} = T;    
p{6} = tspan;
p{7} = A_P;

odefun = @(t,y) slirp_function(t,y,p);

y0(1) = B_I; 
y0(2) = P_I; 
y0(3) = S_I; 
y0(4) = L_I;
y0(5) = I_I; 
y0(6) = R_I; 

[tspan,y] = rk4(odefun,tspan,y0');

B = y(:,1);
P = y(:,2);
S = y(:,3);
L = y(:,4);
I = y(:,5);
R = y(:,6);

subplot(2, 3, k+3);
plot(tspan,P/A_P,'-g',tspan,B/A_P,'--m',tspan,S,'-.k',tspan,L,'--c',tspan,I,':b',tspan,R,'-.r','LineWidth',2);
title('mu_L_min = ', mu_L_min);
xlabel('days (t)');
ylabel('Population (fraction of initial)')
legend('total pop','susceptible pop','berry pop','latent pop ','infected pop','removed pop','Location','northwest')
xlim([0,61]);
hold on

end

function T_B = solve(T)
    if (T<=0)
        T_B = 0;
    elseif (T<35)
        T_B = 0.000241*T^2.06737*(35 - T)^0.72859;
    else
        T_B = 0;
    end
end

function [dydt] = slirp_function(t,y,p)
    B_max = p{1};
    mu_L = p{2};
    mu_I = 1/p{3}; 
    e = p{4};
    T = p{5};
    day = p{6};
    A = p{7};
    B = y(1);
    P = y(2);
    S = y(3);
    L = y(4);
    I = y(5);
    R = y(6);

    if(ceil(t)==floor(t))
        T_new = T(t);
        day_new = day(t);
        mul_new = mu_L(t);
    else 
        t = floor(t);
        T_new = 0.5*(T(t)+T(t+1));
        day_new = 0.5*(day(t)+day(t+1));
        mul_new = 0.5*(mu_L(t)+mu_L(t+1));
    end
    beta = B_max*solve(T_new); 
    TE = -0.35968 +0.10789*T_new-0.00214*T_new^2;
    dP_Bdt = (0.1724*B-0.0000212*B*B)*TE;
    dP_dt = (1.33*(day_new+30))*TE+dP_Bdt; 
    dSdt = -beta*S*I+dP_dt/A; 
    dLdt = beta*S*I-mul_new*L+e;  
    dIdt = mul_new*L-mu_I*I;      
    dRdt = mu_I*I;          
    dydt = [dP_Bdt; dP_dt; dSdt; dLdt; dIdt; dRdt];
    dydt = dydt';

end

function [t, y] = rk4(odefun, tspan, y0)
    num_steps = length(tspan);
    num_states = length(y0);
    t = tspan;
    y = zeros(num_steps, num_states);
    y(1,:) = y0(:);
for i = 1:(num_steps-1)
         h = tspan(i+1) - tspan(i);
         k1 = odefun((i), y(i,:));
         k2 = odefun((i)+0.5, y(i,:) + 0.5*k1*h);
         k3 = odefun((i)+0.5, y(i,:) + 0.5*k2*h);
         k4 = odefun((i+1), y(i,:) + k3*h);
         slope = (k1+ 2*k2 + 2*k3 + k4) / 6;
         y(i+1,:) = y(i,:) + h*slope;
         t(i+1) = tspan(i+1);
end 
end 


% Using the nine plots you created in e) and f) discuss the impact of 𝛽^>_ and 𝜇H,^ij
%on the development of an epidemic. Comment if your conclusion on the impact of 𝛽^>_
%is different than what was observed with a constant 𝛽 in Lab 8. Also comment on how the
%introduction of a more realistic plant growth model and the inclusion of real air temperature
%data impact the development of the epidemic?

% Based on the 7 plots that are shown we can distinctively show the
% differences with different Beta and mu_L_min values on the development of
% an epidemic

% for a change in Beta, we can see a similar peak in the berry population
% but different levels of decline. The infected population seems to end up
% at the same point at t = 61 but has different slopes as it fluctuates and
% increases, this same point applies to the latent population, the
% increasing slop fluctuates around t = 40 when beta increases

% For the change in the min mu_L we can see consistency in the total
% population as well as the susceptible population, the values that
% fluctuate are once again the berry population which over an increase in
% mu_L has a greater peak population. The infected population, removed population, and the latent population decreases
% over time.

% the difference of the beta value from Lab 08 was this graph provides a
% higher population fraction and shows the model on a population increase
% while Lab 08 showed a peak and decline


% More realistic plant growth fluctuates during the different seasons of
% the year as well as the temperature that cannot necessarily be predicted
% through this model. The epidemic with this seperate model would change
% the fluctuation of the numbers

