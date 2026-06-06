%Load Experimental Voltage Traces at both pacing frequencies 
data = readmatrix('voltage_trace.csv') ;
data2 = readmatrix('voltage_trace2.csv') ;
t_experimental = data(:, 1) ;
V_experimental = data(:, 2) ; 
t_experimental2 = data2(:, 1) ;
V_experimental2 = data2(:, 2) ;

%Initial Conditions for HH Model
t_span = [0, 100] ; %Run for 100 seconds 
y0 = [-65, 0.05, 0.6, 0.32] ; %V_0, m_0, h_0, n_0 

%Set bounds and options
lb = [60, 18, 0.15] ; %gNa, gK, gL lower bounds
ub = [180, 54, 0.45] ; %gNa, gK, gL upper bounds
options = optimoptions('ga', 'PopulationSize',200, 'MaxGenerations', 200, ...
    'CrossoverFraction', 0.5, 'Display', 'iter', PlotFcn='gaplotbestf') ;
%Function for extra arguments
fitnessfcn = @(p) evaluate(p, t_span, y0, ...
t_experimental, V_experimental, t_experimental2, V_experimental2);
%Call GA 
[best_params, best_MSE] = ga(fitnessfcn, 3, ...
    [], [], [], [], lb, ub, [], options); 
%Display results
gNa_found = best_params(1) ;
gK_found = best_params(2) ;
gL_found = best_params(3) ; 

fprintf('\nTrue:  gNa=120, gK=36, gL=0.3\n');
fprintf('Found: gNa=%.2f, gK=%.2f, gL=%.4f\n', ...
        gNa_found, gK_found, gL_found);
fprintf('Final MSE: %.4f\n', best_MSE);

%Plot recovered vs. experimental trace
pars_best = [10, gNa_found, gK_found, gL_found, 50, -77, -54.4, 1.0] ;
[t_best, y_best] = ode15s(@derivatives, t_span, y0, [], pars_best) ; 

figure
hold on
plot(t_experimental, V_experimental, '-b') 
plot(t_best, y_best(:,1), '-r') 
xlabel('Time (s)') 
ylabel('Voltage (mV)')
legend('Experimental', 'Simulated') 

%MSE Function for Simulated Protocols 
function MSE = evaluate(p, t_span, y0, ...
    t_experimental, V_experimental, t_experimental2, V_experimental2)
    %Unpack params that vary 
    gNa = p(1) ;
    gK = p(2) ; 
    gL = p(3) ;
    %Protocol 1 (I_app 10)
    pars1 = [10, gNa, gK, gL, 50, -77, -54.4, 1.0] ;
    [t1, y1] = ode15s(@derivatives, t_span, y0, [], pars1) ;
    V_sim1 = interp1(t1, y1(:,1), t_experimental);
    MSE1 = mean((V_sim1 - V_experimental).^2) ;
    %Protocol 2 (I_app 15)
    pars2 = [15, gNa, gK, gL, 50, -77, -54.4, 1.0] ;
    [t2, y2] = ode15s(@derivatives, t_span, y0, [], pars2) ;
    V_sim2 = interp1(t2, y2(:,1), t_experimental2);
    MSE2 = mean((V_sim2 - V_experimental2).^2) ;
    %Return MSE 
    MSE = MSE1 + MSE2 ;
end

%Derivatives Function for HH 
function output = derivatives(t, y, p) 
    derivs = zeros(size(y)) ;
    %Unpack parameters 
    I_app = p(1) ;
    gNa = p(2) ;
    gK = p(3) ;
    gL = p(4) ; 
    ENa = p(5) ; 
    EK = p(6) ;
    EL = p(7) ;
    Cm = p(8) ;
    %Unpack state variables
    V = y(1) ;
    m = y(2) ; 
    h = y(3) ; 
    n = y(4) ;
    %calculations for rate functions
    alpha_m_V = 0.1 * (V + 40) / (1 - exp(-(V + 40)/10)) ;
    beta_m_V = 4.0 * exp(-(V+65)/18) ;
    alpha_h_V = 0.07 * exp(-(V+65)/20) ;
    beta_h_V = 1 / (1 + exp(-(V+35)/10)) ;
    alpha_n_V = 0.01 * (V + 55) / (1 - exp(-(V+55)/10)) ;
    beta_n_V = 0.125 * exp(-(V + 65)/80) ;
    %calculations for currents
    INa = gNa * m.^3 * h * (V - ENa) ;
    IK = gK * n.^4 * (V - EK) ;
    IL = gL * (V - EL) ;
    %compute derivatives
    derivs(1) = (1/Cm) * (I_app - INa - IK - IL) ;
    derivs(2) = alpha_m_V * (1 - m) - beta_m_V * m ;
    derivs(3) = alpha_h_V * (1 - h) - beta_h_V * h ;
    derivs(4) = alpha_n_V * (1 - n) - beta_n_V * n ;
    output = derivs ;
end