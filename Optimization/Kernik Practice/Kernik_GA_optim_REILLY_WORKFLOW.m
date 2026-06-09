load ICs_baseline
load baseline_parameter_inputs

%Parallel computing 
if isempty(gcp('nocreate'))
    parpool('local') ;
end

%Load Experimental Voltage Traces generated from Kernik Model at both pacing frequencies 
%CL of 1000 and 2000
data = readmatrix('Kernik_values_1000.csv') ;
data2 = readmatrix('Kernik_values_2000.csv') ;
t_experimental_1 = data(:, 1) ;
Ca_experimental_1 = data(:, 3) ; 
t_experimental_2 = data2(:, 1) ;
Ca_experimental_2 = data2(:, 3) ; 

%Initial Conditions for Kernik
t_span = [0, 10.6e4] ; %Run for 106 seconds
%Set bounds and options
lb = exp(-1) * ones(1,16) ;
ub = exp(1)  * ones(1,16) ;
options = optimoptions('ga', ...
    'PopulationSize', 300, ...   
    'MaxGenerations', 50, ...    
    'EliteCount', 5, ...
    'CrossoverFraction', 0.8, ...
    'UseParallel', true, ...
    'Display', 'iter') ;

%Function for extra arguments
fitnessfcn = @(p) evaluate(p, baseline_parameter_inputs, Y_init, t_span,...
t_experimental_1, Ca_experimental_1, ...
t_experimental_2, Ca_experimental_2) ;

%Call GA 
[best_params, best_MSE] = ga(fitnessfcn, 16, ...
    [], [], [], [], lb, ub, [], options); 

%Display Results
param_names = {'IK1','IKr','IKs','Ito','ICaL','ICaT','INa',...
               'If','NCX','SERCA','RyR','Jleak','INaK',...
               'IbNa','IbCa','IpCa'} ;

true_scales = ones(1,16) ;
true_scales(2)  = 0.7 ;
true_scales(5)  = 1.2 ;
true_scales(8)  = 1.3 ;
true_scales(9)  = 0.5 ;
true_scales(10) = 0.5 ;
true_scales(11) = 1.5 ;
true_scales(15) = 1.3 ;
true_scales(16) = 0.6 ;

fprintf('\nParameter Recovery:\n') ;
for i = 1:16
    fprintf('%s: True=%.2f Found=%.4f\n', ...
            param_names{i}, true_scales(i), best_params(i)) ;
end

%Simulation
pars_best = baseline_parameter_inputs ;
pars_best(1:16) = best_params ;  

options_ODE = odeset('MaxStep',1,'InitialStep',2e-2) ;
[t_best, Y_best] = ode15s(@ipsc_function, t_span, ...
                          Y_init, options_ODE, pars_best) ;

%Plot Results
figure
hold on
last_6_best = t_best >= 100000 ;
plot(t_experimental_1, Ca_experimental_1, '-b', 'LineWidth', 2)
plot(t_best(last_6_best), Y_best(last_6_best,3), '--r', 'LineWidth', 2)
xlabel('Time (ms)') ; ylabel('Calcium (mM)')
legend('Experimental', 'Recovered')

%SAD Function for Simulated Protocols 
function SAD = evaluate(scales, base_params, Y_init, t_span, ...
    t_exp_1, Ca_exp_1, t_exp_2, Ca_exp_2)
    %Unpack scales
    pars = base_params ;
    pars(1:16) = scales ;
    %Protocol 1 - CL of 800
    options_ODE = odeset('MaxStep',1,'InitialStep',2e-2) ;
    try
        [t1, y1] = ode15s(@ipsc_function, t_span, Y_init, options_ODE, pars) ;
        %Extract last 6 seconds
        last_6_idx = t1 >= 10e4 ;
        t1_last  = t1(last_6_idx) ;
        Ca_sim_1 = y1(last_6_idx, 3) ;
        V_sim_1  = y1(last_6_idx, 1) ;

        %Biomarkers 
        AP_amp   = max(V_sim_1) - min(V_sim_1) ;
        RMP      = min(V_sim_1) ;
        dVdtmax  = max(diff(V_sim_1)./diff(t1_last)) ;
        Ca_amp   = max(Ca_sim_1) - min(Ca_sim_1) ;
        Ca_diast = min(Ca_sim_1) ;

        % Penalties
        P_AP   = compute_penalty(AP_amp,   90,     125,    105,  600) ;
        P_RMP  = compute_penalty(RMP,      -80,    -55,    -75,  600) ;
        P_dVdt = compute_penalty(dVdtmax,  10,     150,    25,   600) ;
        P_Ca   = compute_penalty(Ca_amp,   1.5e-4, 1e-3,   5e-4, 600) ;
        P_Cd   = compute_penalty(Ca_diast, 1e-4,   3.5e-4, 2e-4, 600) ;
        P_total = P_AP + P_RMP + P_dVdt + P_Ca + P_Cd ;

        % Normalize simulation
        Ca_sim_norm = (Ca_sim_1 - min(Ca_sim_1)) / (max(Ca_sim_1) - min(Ca_sim_1)) ;
        % Normalize experimental
        Ca_exp_norm = (Ca_exp_1 - min(Ca_exp_1)) / (max(Ca_exp_1) - min(Ca_exp_1)) ;
        Ca_interp = interp1(t1_last, Ca_sim_norm, t_exp_1) ;
        SAD_1 = sum(abs(Ca_interp - Ca_exp_norm)) ;
    catch
        SAD = 999999 ; return ;
    end
    %Protocol 2 - CL of 500
    try
        options_ODE = odeset('MaxStep',1,'InitialStep',2e-2) ;
        [t2, y2] = ode15s(@ipsc_function_2000, t_span, Y_init, options_ODE, pars) ;
        last_6_idx = t2 >= 10e4  ;
        t2_last  = t2(last_6_idx) ;
        V_sim_2  = y2(last_6_idx, 1) ;
        Ca_sim_2 = y2(last_6_idx, 3) ;

        %Biomarkers 2
        AP_amp2   = max(V_sim_2) - min(V_sim_2) ;
        RMP2      = min(V_sim_2) ;
        dVdtmax2  = max(diff(V_sim_2)./diff(t2_last)) ;
        Ca_amp2   = max(Ca_sim_2) - min(Ca_sim_2) ;
        Ca_diast2 = min(Ca_sim_2) ;

        % Penalties
        P_AP_2   = compute_penalty(AP_amp2,   90,     125,    105,  600) ;
        P_RMP_2  = compute_penalty(RMP2,      -80,    -55,    -75,  600) ;
        P_dVdt_2 = compute_penalty(dVdtmax2,  10,     150,    25,   600) ;
        P_Ca_2   = compute_penalty(Ca_amp2,   1.5e-4, 1e-3,   5e-4, 600) ;
        P_Cd_2   = compute_penalty(Ca_diast2, 1e-4,   3.5e-4, 2e-4, 600) ;
        P_total_2 = P_AP_2 + P_RMP_2 + P_dVdt_2 + P_Ca_2 + P_Cd_2 ;
        % Normalize simulation
        Ca_sim_norm_2 = (Ca_sim_2 - min(Ca_sim_2)) / (max(Ca_sim_2) - min(Ca_sim_2)) ;
        % Normalize experimental
        Ca_exp_norm_2 = (Ca_exp_2 - min(Ca_exp_2)) / (max(Ca_exp_2) - min(Ca_exp_2)) ;
        Ca_interp_2 = interp1(t2_last, Ca_sim_norm_2, t_exp_2) ;
        SAD_2 = sum(abs(Ca_interp_2 - Ca_exp_norm_2)) ;
    catch
        SAD = 999999 ; return ;
    end

    lambda = 0.1 ;
    lasso_penalty = lambda * sum(abs(scales - 1)) ;

    %Return SAD 
    SAD = SAD_1 + SAD_2 + P_total + P_total_2 + lasso_penalty ;
end

%penalty function
function P = compute_penalty(value, Min, Max, Avg, gamma)
    if value > Min && value < Max
        P = 0 ;
    else
        P = gamma * min(abs((value-Min)/Avg), ...
                        abs((value-Max)/Avg)) ;
    end
end