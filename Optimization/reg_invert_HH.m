%Load Experimental Voltage Traces at both pacing frequencies 
data = readmatrix('voltage_trace.csv') ;
data2 = readmatrix('voltage_trace2.csv') ;
t_experimental = data(:, 1) ;
V_experimental = data(:, 2) ; 
t_experimental2 = data2(:, 1) ;
V_experimental2 = data2(:, 2) ;

%Using Linear Regression for Inversion instead of GA 
n_sims = 1000 ; 

%Generate random params for currents 
lb = [60, 18, 0.15] ; %gNa, gK, gL lower bounds
ub = [180, 54, 0.45] ; %gNa, gK, gL upper bounds
gNa_array = lb(1) + rand(n_sims, 1) * (ub(1) - lb(1)) ;
gK_array  = lb(2) + rand(n_sims, 1) * (ub(2) - lb(2)) ;
gL_array  = lb(3) + rand(n_sims, 1) * (ub(3) - lb(3)) ;

%Prellaocate Storage 
biomarkers_matrix = zeros(n_sims, 10) ;
pars_matrix = zeros(n_sims, 3) ; 

%Initial Conditions for HH Model
t_span = [0, 100] ; %Run for 100 seconds 
y0 = [-65, 0.05, 0.6, 0.32] ; %V_0, m_0, h_0, n_0 

%Run 1000 sims with varied parameters 
for i = 1:n_sims 
    %Protocol 1 
    pars = [10, gNa_array(i), gK_array(i), gL_array(i), 50, -77, -54.4, 1.0] ;
    [t,y] = ode15s(@derivatives, t_span, y0, [], pars) ;
    V = y(:,1) ;
    %Protocol 2
    pars2 = [15, gNa_array(i), gK_array(i), gL_array(i), 50, -77, -54.4, 1.0] ;
    [t2, y2] = ode15s(@derivatives, t_span, y0, [], pars2) ;
    V2 = y2(:,1) ;
    %Extract biomarkers 
    biomarkers_matrix(i, 1) = max(V) ; %Peak Voltage
    biomarkers_matrix(i, 2) = min(V) ; %Min Voltage
    biomarkers_matrix(i, 3) = max(diff(V)./diff(t)) ; %Max upstroke
    biomarkers_matrix(i, 4) = V(end) ; %RMP
    biomarkers_matrix(i, 5) = compute_APD(V, t, 0.1) ; %APD90 
    biomarkers_matrix(i, 6)  = max(V2) ;
    biomarkers_matrix(i, 7)  = min(V2) ;
    biomarkers_matrix(i, 8)  = max(diff(V2)./diff(t2)) ;
    biomarkers_matrix(i, 9)  = V2(end) ;
    biomarkers_matrix(i, 10) = compute_APD(V2, t2, 0.1) ;
    %Extract params 
    pars_matrix(i,:) = [gNa_array(i), gK_array(i), gL_array(i)] ;
end

%NORMALIZE FOR REGRESSION
bio_mean = mean(biomarkers_matrix) ;
bio_std = std(biomarkers_matrix) ;
par_mean = mean(pars_matrix) ;
par_std = std(pars_matrix) ;

%Normalize matrices
biomarkers_norm = (biomarkers_matrix - bio_mean) ./ bio_std ;
pars_norm = (pars_matrix - par_mean) ./ par_std ;

%Fit Regression
[XL, YL, XS, YS, BETA, PCTVAR] = plsregress(biomarkers_norm, pars_norm, 3) ;
fprintf('PLS explains %.1f%% of parameter variance\n', ...
        sum(PCTVAR(2,:)) * 100) ;

% Extract biomarkers from experimental trace
exp_bio = [max(V_experimental), ...
           min(V_experimental), ...
           max(diff(V_experimental)./diff(t_experimental)), ...
           V_experimental(end), ...
           compute_APD(V_experimental, t_experimental, 0.1), ...
           max(V_experimental2), ...
           min(V_experimental2), ...
           max(diff(V_experimental2)./diff(t_experimental2)), ...
           V_experimental2(end), ...
           compute_APD(V_experimental2, t_experimental2, 0.1)] ;
% Normalize using TRAINING statistics
exp_bio_norm = (exp_bio - bio_mean) ./ bio_std ;
% Predict parameters
pars_norm_pred = [1, exp_bio_norm] * BETA ;
% Unnormalize
pars_recovered = pars_norm_pred .* par_std + par_mean ;

% Display
fprintf('\nTrue:  gNa=120, gK=36, gL=0.3\n') ;
fprintf('Found: gNa=%.2f, gK=%.2f, gL=%.4f\n', ...
        pars_recovered(1), pars_recovered(2), pars_recovered(3)) ;

%FUNCTIONS
%CALCULATE APD 
function APD = compute_APD(voltage, time, repol_percent) 
    V_thresh = voltage(1) + repol_percent * (max(voltage) - voltage(1)) ; %Resting voltage + appropriate % of max - resting voltage
    %Calculate APD based on crossing threshold up or down 
    above_thresh = voltage > V_thresh ;
    crossings = diff(above_thresh) ;
    up_crossings = find(crossings == 1) ; %upstroke
    down_crossings = find(crossings == -1) ; %repol
    %Check if no APs fire 
    if isempty(up_crossings) || isempty(down_crossings)
        APD = 0 ;
        return ; 
    end
    n_pairs = min(length(up_crossings), length(down_crossings)) ;
    AP_diff = time(down_crossings(1:n_pairs)) - time(up_crossings(1:n_pairs)) ;
    APD = mean(AP_diff) ; 
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