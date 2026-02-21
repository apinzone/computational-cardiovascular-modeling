%expdecay model with ode15s

%parameters
g = 9.8;
k = 0.5;
params = [g , k];
%initial conditions
t_span = [0 ; 10];
y0 = [0 ; 100]; % v, y
options = [];

%solvercall
[t, y] = ode15s(@derivatives, t_span, y0, options, params);

%plot 
plot(t, y(:,1));   %velocity over time

%equations
function output = derivatives(t, y, p)
    %unpack parameters
    g = p(1);
    k = p(2);
    %state variables
    v = y(1);

    %compute derivatives
    dvdt = -g - k * v;
    dydt = v;
    %return array of derivatives 
    output = [dvdt ; dydt];
end

