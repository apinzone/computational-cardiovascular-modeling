%Let's get back on track , Fitzhugh-Nagumo


%Initial Conditions
t_span = [0, 200] ; %Run for 200s
y0 = [-1, -0.5]; %V0, w0 

%Parameters
I_range = -0.5:0.3:1.5 ; %Vary I from -0.5 to 1.5 by increments of 0.3
bifurcation_array_V_max = zeros(size(I_range)) ;
bifurcation_array_V_min = zeros(size(I_range)) ;
%Solver call 
for i = 1:length(I_range)
    pars = [0.7, 0.8, 0.08, I_range(i)] ; %a, b, epsilon, I
    [t, y] = ode15s(@derivatives,t_span, y0, [], pars) ;
    last_50_index = round(0.5 * length(t)) ;
    bifurcation_array_V_max(i) = max(y(last_50_index:end,1)); 
    bifurcation_array_V_min(i) = min(y(last_50_index:end,1));
end

%Bifurcation diagram with max and min voltages across range of I values
figure
hold on
plot(I_range, bifurcation_array_V_max, '-b') ;
xlabel('I Values')
plot(I_range, bifurcation_array_V_min, '-r') ;
legend('V_max', 'V_min')

%Voltage over time plot
figure
hold on
plot(t, y(:,1));
xlabel('Time')
ylabel('Voltage (mV)')

%Derivatives Function
function output = derivatives(t, y, p) 
    derivs = zeros(size(y)) ;
    %Unpack State Variables
    V = y(1) ;
    w = y(2) ;

    %Unpack params 
    a = p(1) ;
    b = p(2) ;
    epsilon = p(3) ;
    I = p(4) ; 

    %Compute derivatives
    derivs(1) = V - V.^3/3 - w + I ; %dVdt
    derivs(2) = epsilon * (V + a - b * w) ; %dwDt
    output = derivs ;
end


