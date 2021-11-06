function [battery_required, Max_flight_time] = Energy_consumption_modelization(t,v)
%Output 
%battery_tr = percentage of battery required to complete the segments given
%the time t and the velocity v

%Input 
% t = time of each segments (length(idxs) * 1)
% v =  velocity horizontal o verticale

% v can be referred to the orizontal or vertical velocity:
% if horizontal = 1  >>> evaluate the energy used in the orizontal flight
% if horizontal = 0  >>> evaluate the energy used in the vertical flight
% Battery Capacity and Energy Consumption evaluation for mavic pro 2 enterprsise

horizontal = 1; % evaluate orizontal flight 
% Dynamic related parameters 
m = 734 * 10^(-3); % [Kg]
g = 9.89; % [m/s^2] gravity acceleration
A = (214* 10^(-3)) * (91 * 10^(-3)) * (82 * 10^(-3)); % [m^3] horizontal cross section area
Cd = 1; % drag coefficient
rho = 1.29; % [Kg/m^3] air density
Ct = 0.03; % thrust coefficent 
r = 25 * 10^(-2); % [m] radius of propellers 
N = 4; % number of propellers
Ap = pi*r^2; % disk area of propellers

% v = 8.3; %[m/s] desired horizontal velocity
% t = 1200; %[sec] fligth time in seconds

% Energy related parameters
cut_off_p = 30; % cut off percentage
Cmax = 59.29 * 60 * 60; % [W * s] maximum capacity of mavic 2 battery
% Cmin = cut_off_p * Cmax / 100;
% C = Cmax - Cmin;
V = 15.4; % battery voltage 


% Evaluation of Angular velocity of the propellers for a horizontal flight
% with constant velocity v = 8.3 m/s
P = zeros(length(t),1);
E = zeros(length(t),1);
battery_required = zeros(length(t),1); % initialize the battery vector (same size of the time vector) that will contains the battery required for each link 

% if we are dealing with the orizontal flight evaluation 
if horizontal == 1
    
    n1= (((m * g)^2 + (1/2 * rho * A * Cd * v^2)^2)^(1/4));
    n2 = (r^2 * N * rho * A * Ct)^(1/2);
    w = n1/n2;
% if we are dealing with the vertical flight evaluation    
else
    % calcola energia spesa in volo verticale
    n1= ((2*(m * g) + (rho * A * Cd * v^2))^(1/2));
    n2 = (r^2 * N * rho * A * Ct)^(1/2);
    w = n1/n2;
end

% Power Consumption curve given the angular velocity with parameters evaluated from experimental results
P = 2.258*10^(-7) * w^3 + 3.866 * 10^(-5) * w^2 + 5.137 * 10^(-3) * w + 2.616; %[Watt]

for ii = 1 : length(t) % fill the battery vector
% Energy consumption [Watt/second]
E(ii) = P * t(ii);

Max_flight_time = Cmax/P; % At a given velocity v

battery_required(ii) = E(ii) * 100 / Cmax;

% Normalization between zero and 1
battery_required(ii) = battery_required(ii) / 100;

end
end
% %Evaluation of the power consumption in function of velocity
% v_test = [0 : 0.01 : 90]; %m/s
% for i = 1 : length(v_test)
%     n1= (((m * g)^2 + ((1/2) * rho * A * Cd * v_test(i)^2)^2)^(1/4));
%     n2 = (r^2 * N * rho * A * Ct)^(1/2);
%     w_test(i) = n1/n2;
%     P(i) = 2.258*10^(-7) * w_test(i)^3 + 3.866 * 10^(-5) * w_test(i)^2 + 5.137 * 10^(-3) * w_test(i) + 2.616; 
% end
% 
% 
% figure 
% plot(P)
% xlim([0 90]);ylim([90 100])


%xlim([0 10+2]);ylim([0 10+2])
