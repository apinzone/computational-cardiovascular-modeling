


%parameters
alpha = 1; %prey birth rate
beta = 0.1; %predation rate
delta = 0.075; %predator growth rate
gamma = 1.5; %predator death rate

%initial conditions
P = 10; %starting prey pop
Q = 5.0; %starting predator population
t = 0; 
dt = 0.01;
t_end = 50;

%array creation
prey = [P];
predator = [Q];
time_array = [t];

%solver loop
while t<t_end
    %compute derivatives
    dp = dpdt(P, alpha, beta, Q) * dt;
    dq = dqdt(P, delta, gamma, Q) * dt;
    %update states
    P = P + dp;
    Q = Q + dq;
    t = t + dt;
    %append arrays
    prey(end + 1) = P;
    predator(end + 1) = Q;
    time_array(end + 1) = t; % append current time to time array
end

plot(time_array, prey)
hold on 
plot(time_array, predator)
xlabel('Time')
ylabel('Population')

%define equations
function dp = dpdt(P, alpha, beta, Q)
    dp = alpha * P - beta * P * Q;
end 

function dq = dqdt (P, delta, gamma, Q)
    dq = delta * P * Q - gamma * Q;
end