
%parameters
I = 0.5; %external stimulus current

%initial conditions
t = 0;
t_end = 100;
dt = 0.01;
V = -1.0;
W = 0.0;

%arrays
voltage_array = [V];
w_array = [W];
time_array = [t];

%solver loop
while t<t_end
    %compute derivatives
    dV = dVdt(V, W, I) * dt;
    dW = dWdt(V, W) * dt;
    %Update state variables
    V = V + dV;
    W = W + dW;
    t = t + dt;
    %Store updated values in arrays
    voltage_array = [voltage_array; V];
    w_array = [w_array; W];
    time_array = [time_array; t];
end 

%Plot results
plot(time_array, voltage_array);
xlabel('Time')
ylabel('Voltage')

figure
plot(w_array, voltage_array);


%Define Equations
function dV = dVdt(V, W, I)
    dV = V - (V.^3/3) - W + I;
end

function dW = dWdt(V, W)
    dW = (V + 0.7 - 0.8 * W)/12.5;
end






