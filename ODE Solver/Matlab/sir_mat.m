%implement SIR model


%parameters [mu, beta, gamma]
pars = [0.02; 0.5; 0.1];

%initial conditions [S, I, R]
y0 = [990; 10; 0];
t_span = [0; 200];

%call solver 
[t, y] = ode15s(@derivatives, t_span, y0, [], pars);
plot(t, y(:,1)); %PLOT S 
hold on;
plot(t, y(:,2)); %PLOT I
hold on;
plot(t, y(:,3)); %PLOT R
hold on;
legend('S', 'I', 'R');
xlabel('Time (s)')
ylabel('Population Parameters')



%derivatives function
function output = derivatives(t, y, p)
    derivs = zeros(size(y));
    %unpack params
    mu = p(1);
    beta = p(2);
    gamma = p(3);

    %unpack state variables
    S = y(1);
    I = y(2);
    R = y(3);
    %calculate N
    N = S + I + R;

    %compute derivatives
    derivs(1) = mu * N - beta * (S * I)/ N - (mu * S);
    derivs(2) = beta * (S * I) / N - (gamma * I) - (mu * I);
    derivs(3) = (gamma * I) - (mu * R);

    output = derivs;
end 


