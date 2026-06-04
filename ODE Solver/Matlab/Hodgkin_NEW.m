%Parameters and Initial Conditions
y0 = [-65, 0.05, 0.6, 0.32] ; %V_0, m_0, h_0, n_0 
t_span = [0, 100] ; %Run for 100s. 
I_app_range = linspace(0, 20, 50) ;
%Prellocate arrays for bifurc diagrams
V_max_array = zeros(size(I_app_range)) ; 
V_min_array = zeros(size(I_app_range)) ; 

%Solver Call With Parameter Sweep
for i = 1:length(I_app_range) 
    %i_app, gNa, gK, gL, ENa, EK, EL, Cm
    pars = [I_app_range(i), 120, 36, 0.3, 50, -77, -54.4, 1.0] ; 
    [t, y] = ode15s(@derivatives, t_span, y0, [], pars) ;
    t_idx_settled = round(0.5 * length(t)) ;
    V_max_array(i) = max(y(t_idx_settled:end,1)) ; 
    V_min_array(i) = min(y(t_idx_settled:end,1)) ; 
end

%Bifurcation Plot 
figure 
hold on
plot(I_app_range, V_min_array, '-b') ;
plot(I_app_range, V_max_array, '-r') ;
xlabel('I_app')
ylabel('Voltage (mV)')
legend('V min', 'V max') ; 

%Voltage over time
figure
plot(t, y(:,1)) ; 

%Derivatives Function
function output = derivatives(t, y, p) 
    derivs = zeros(size(y)) ;
    %unpack parameters
    I_app = p(1) ;
    gNa = p(2);
    gK = p(3) ;
    gL = p(4) ;
    ENa = p(5) ;
    EK = p(6) ;
    EL = p(7) ;
    Cm = p(8) ; 

    %unpack state variables 
    V = y(1) ;
    m = y(2) ;
    h = y(3) ;
    n = y(4) ;

    %Calculations for Rate Functions
    alpha_m_V = 0.1 * (V + 40) / (1 - exp(-(V + 40)/10)) ;
    beta_m_V = 4.0 * exp(-(V+65)/18) ;
    alpha_h_V = 0.07 * exp(-(V+65)/20) ;
    beta_h_V = 1 / (1 + exp(-(V+35)/10)) ;
    alpha_n_V = 0.01 * (V + 55) / (1 - exp(-(V+55)/10)) ;
    beta_n_V = 0.125 * exp(-(V + 65)/80) ; 

    %Calculations for Currents 
    INa = gNa * m.^3 * h * (V - ENa) ;
    IK = gK * n.^4 * (V - EK) ;
    IL = gL * (V - EL) ;

    %compute derivatives and append deriv array
    derivs(1) = (1/Cm) * (I_app - INa - IK - IL) ; %dV/dt
    derivs(2) = alpha_m_V * (1 - m) - beta_m_V * m ;%dm/dt
    derivs(3) = alpha_h_V * (1 - h) - beta_h_V * h ; %dh/dt
    derivs(4) = alpha_n_V * (1 - n) - beta_n_V * n ; %dn/dt

    %Return derivatives
    output = derivs ; 
end
