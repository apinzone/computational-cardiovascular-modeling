%implement pendulum model with damper 

%define params [g, L, b]
pars = [9.8; 1.0; 0.5];

%set initial conditions
t_span = [0;20]; % simulate for 20 seconds
y0 = [pi/2, 0]; %theta, omega

%call solver 
[t, y] = ode15s(@derivatives, t_span, y0, [], pars);

theta_plot = y(:, 1);
omega_plot = y(:, 2);
%plot
plot(t, theta_plot)
hold on;
plot(t, omega_plot)

xlabel('time (s)')

%derivatives function
function output = derivatives(t, y, p)
    %blank array for derivs
    derivatives = zeros(size(y));
    %unpack paramaeters
    g = p(1);
    L = p(2);
    b = p(3);
    %state variables
    theta = y(1);
    omega = y(2);

    %compute derivatives but update simultaneously
    derivatives(1) = omega; %theta
    derivatives(2) = -(g/L) * sin(theta) - b * omega; %omega

    %return derivatives
    output = derivatives;
end
