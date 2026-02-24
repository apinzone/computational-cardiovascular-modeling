%parameters [a, b, c, d (burst 1), a, b, c, d (burst 0)]
pars = [0.02; 0.2; -50; 2; 0.02; 0.2; -65; 8];

%initial conditions
y0 = [-65 ; -13]; %V, U
t_span = [0 ; 300];

%call solver 
[t, y] = ode15s(@derivatives, t_span, y0, [], pars, 1);

%plot
plot(t,y(:,1))
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Membrane Potential Over Time');
grid on;

%derivatives function 
function output = derivatives(t, y, p, bm)
    derivs = zeros(size(y));
    %parameters based on condition
    I = 10; %constant
    if bm == 1
        a = p(1);
        b = p(2);
        c = p(3);
        d = p(4);

    else
        a = p(5);
        b = p(6);
        c = p(7);
        d = p(8);

    end 
    %state variables
    v = y(1);
    u = y(2);

    %compute derivatives
    derivs(1) = (0.04*v.^2) + (5 * v) + 140 - u + I;
    derivs(2) = a * ((b*v) - u);

    output = derivs;
end
 