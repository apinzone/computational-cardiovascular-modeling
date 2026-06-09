%Compute APD and CaT from Voltage and Calcium Traces, and Extract Final
%Beat

data = readmatrix('Kernik_values_1000') ;
data2 = readmatrix('Kernik_values_2000') ;
V_1000 = data(:,2) ;
t_1000 = data(:, 1) ;
Cai_1000 = data(:, 3) ;

V_2000 = data2(:,2) ;
t_2000 = data2(:, 1) ;
Cai_2000 = data2(:, 3) ;
t_end = 10.6e4 ; %Simulation runs for 106 s

AP_amp_1000 = max(V_1000) - min(V_1000) ;
AP_amp_2000 = max(V_2000) - min(V_2000) ;
t_eval = linspace(0, 10.6e4, length(V_1000)) ;
dVdtmax_1000 = max(diff(V_1000)./diff(t_1000)) ;
dVdtmax_2000 = max(diff(V_2000)./diff(t_2000)) ;

%APD Metrics - 1000 BCL
[APD_30_1000, num_APs_norm_30_1000] = compute_APD(V_1000, t_1000, 30) ;
[APD_60_1000, num_APs_norm_60_1000] = compute_APD(V_1000, t_1000, 60) ;
[APD_90_1000, num_APs_norm_90_1000] = compute_APD(V_1000, t_1000, 90) ;
fprintf('APD at 30 (1000 BCL) is %.2f\n', APD_30_1000) ;
fprintf('APD at 60 (1000 BCL)is %.2f\n', APD_60_1000) ;
fprintf('APD at 90 (1000 BCL) is %.2f\n', APD_90_1000) ;
fprintf('Action Potentials (1000 BCL) (30) (1000 BCL) Captured: %.2f\n', num_APs_norm_30_1000) ;
fprintf('AP amplitude 1000 BCL: %.2f mV\n', AP_amp_1000) ;
fprintf('dVdtmax 1000 BCL: %.2f V/s\n', dVdtmax_1000) ;

%CaT Metrics - 1000 BCL
[Ca_Peak_1000, Ca_Diastolic_1000, Ca_Amplitude_1000, CaD50_1000] = compute_CaT(Cai_1000, t_1000) ; 
fprintf('Peak Calcium (1000 BCL) is %.2f\n', Ca_Peak_1000) ;
fprintf('Diastolic Calcium (1000 BCL) is %.2f\n', Ca_Diastolic_1000) ;
fprintf('Calcium Ampltiude (1000 BCL) is %.2f\n', Ca_Amplitude_1000) ;
fprintf('CaD50 is (1000 BCL) %.2f\n', CaD50_1000) ;

%Last Beat Metrics - 1000 BCL
last_beat_time_1000 = t_end - 1000 ; 
last_beat_idx_1000 = t_1000 >= last_beat_time_1000 ; 
V_last_beat_1000 = V_1000(last_beat_idx_1000) ; 
t_last_beat_1000 = t_1000(last_beat_idx_1000) ; 
Ca_last_beat_1000 = Cai_1000(last_beat_idx_1000) ; 
RMP_LB_1000 = V_last_beat_1000(1) ; 

%APD Metrics Last Beat - 1000 BCL
[APD_30_LB_1000, num_APs_LB_30_1000] = compute_APD(V_last_beat_1000, t_last_beat_1000, 30) ;
[APD_60_LB_1000, num_APs_LB_60_1000] = compute_APD(V_last_beat_1000, t_last_beat_1000, 60) ;
[APD_90_LB_1000, num_APs_LB_90_1000] = compute_APD(V_last_beat_1000, t_last_beat_1000, 90) ;
fprintf('APD at 30 (LB) (1000 BCL) is %.2f\n', APD_30_LB_1000) ;
fprintf('APD at 60 (LB) (1000 BCL) is %.2f\n', APD_60_LB_1000) ;
fprintf('APD at 90 (LB) (1000 BCL) is %.2f\n', APD_90_LB_1000) ;
fprintf('Action Potentials Captured (1000 BCL): %.2f\n', num_APs_LB_30_1000) ;
fprintf('RMP (LB) (1000 BCL) is %.2f\n', RMP_LB_1000) ;

%CaT Metrics Last Beat - 1000 BCL
[Ca_Peak_LB_1000, Ca_Diastolic_LB_1000, Ca_Amplitude_LB_1000, CaD50_LB_1000] = compute_CaT(Ca_last_beat_1000, t_last_beat_1000) ; 
fprintf('Peak Calcium (LB) (1000 BCL) is %.2f\n', Ca_Peak_LB_1000) ;
fprintf('Diastolic Calcium (LB) (1000 BCL) is %.2f\n', Ca_Diastolic_LB_1000) ;
fprintf('Calcium Ampltiude (LB) (1000 BCL) is %.2f\n', Ca_Amplitude_LB_1000) ;
fprintf('CaD50 (LB) (1000 BCL) is %.2f\n', CaD50_LB_1000) ;

%APD Metrics - 2000 BCL
[APD_30_2000, num_APs_norm_30_2000] = compute_APD(V_2000, t_2000, 30) ;
[APD_60_2000, num_APs_norm_60_2000] = compute_APD(V_2000, t_2000, 60) ;
[APD_90_2000, num_APs_norm_90_2000] = compute_APD(V_2000, t_2000, 90) ;
fprintf('APD at 30 (2000 BCL) is %.2f\n', APD_30_2000) ;
fprintf('APD at 60 (2000 BCL)is %.2f\n', APD_60_2000) ;
fprintf('APD at 90 (2000 BCL) is %.2f\n', APD_90_2000) ;
fprintf('Action Potentials (30) (2000 BCL) Captured: %.2f\n', num_APs_norm_30_2000) ;
fprintf('AP amplitude 2000 BCL: %.2f mV\n', AP_amp_2000) ;
fprintf('dVdtmax 2000 BCL: %.2f V/s\n', dVdtmax_2000) ;

%CaT Metrics - 2000 BCL
[Ca_Peak_2000, Ca_Diastolic_2000, Ca_Amplitude_2000, CaD50_2000] = compute_CaT(Cai_2000, t_2000) ; 
fprintf('Peak Calcium (2000 BCL)is %.2f\n', Ca_Peak_2000) ;
fprintf('Diastolic Calcium (2000 BCL) is %.2f\n', Ca_Diastolic_2000) ;
fprintf('Calcium Ampltiude (2000 BCL) is %.2f\n', Ca_Amplitude_2000) ;
fprintf('CaD50 (2000 BCL) is %.2f\n', CaD50_2000) ;

%Last Beat Metrics - 2000 BCL
last_beat_time_2000 = t_end - 2000 ; 
last_beat_idx_2000 = t_2000 >= last_beat_time_2000 ; 
V_last_beat_2000 = V_2000(last_beat_idx_2000) ; 
t_last_beat_2000 = t_2000(last_beat_idx_2000) ; 
Ca_last_beat_2000 = Cai_2000(last_beat_idx_2000) ; 
RMP_LB_2000 = V_last_beat_2000(1) ; 

%APD Metrics Last Beat - 2000 BCL
[APD_30_LB_2000, num_APs_LB_30_2000] = compute_APD(V_last_beat_2000, t_last_beat_2000, 30) ;
[APD_60_LB_2000, num_APs_LB_60_2000] = compute_APD(V_last_beat_2000, t_last_beat_2000, 60) ;
[APD_90_LB_2000, num_APs_LB_90_2000] = compute_APD(V_last_beat_2000, t_last_beat_2000, 90) ;
fprintf('APD at 30 (LB) (2000 BCL) is %.2f\n', APD_30_LB_2000) ;
fprintf('APD at 60 (LB) (2000 BCL) is %.2f\n', APD_60_LB_2000) ;
fprintf('APD at 90 (LB) (2000 BCL) is %.2f\n', APD_90_LB_2000) ;
fprintf('Action Potentials Captured (2000 BCL): %.2f\n', num_APs_LB_30_2000) ;
fprintf('RMP (LB) (2000 BCL) is %.2f\n', RMP_LB_2000) ;

%CaT Metrics Last Beat - 2000 BCL
[Ca_Peak_LB_2000, Ca_Diastolic_LB_2000, Ca_Amplitude_LB_2000, CaD50_LB_2000] = compute_CaT(Ca_last_beat_2000, t_last_beat_2000) ; 
fprintf('Peak Calcium (LB) (2000 BCL) is %.2f\n', Ca_Peak_LB_2000) ;
fprintf('Diastolic Calcium (LB) (2000 BCL) is %.2f\n', Ca_Diastolic_LB_2000) ;
fprintf('Calcium Ampltiude (LB) (2000 BCL) is %.2f\n', Ca_Amplitude_LB_2000) ;
fprintf('CaD50 (LB) (2000 BCL) is %.2f\n', CaD50_LB_2000) ;

%APD Function
function [APD, num_APs] = compute_APD(voltage, time, prcnt)
    V_rest = voltage(1) ; %Baseline voltage
    V_peak = max(voltage) ; %Peak Voltage 
    V_thresh = V_rest + (V_peak - V_rest) * (1 - prcnt/100) ; %Compute voltage at threshold based on user input % repolarization
    above_thresh = voltage > V_thresh ; %Stamp indeces where voltage is above threshold
    crossings = diff(above_thresh) ; %convert to 1s above threshold, -1s below threshold
    up_crossings = find(crossings == 1); %Stamp indeces going up
    down_crossings = find(crossings == -1) ;%Stamp indeces going down
    if isempty(up_crossings) || isempty(down_crossings) %Fallback
        APD = 0 ;
        return ;
    end
    n_pairs = min(length(up_crossings), length(down_crossings)) ;
    AP_diff = time(down_crossings(1:n_pairs)) - time(up_crossings(1:n_pairs)) ;
    APD = mean(AP_diff) ;
    num_APs = length(AP_diff) ;
end

%CaT Function
function [peak, diastolic, amplitude, CaD50] = compute_CaT(calcium, time)
    %convert to nM
    calcium = calcium * 1e6 ;
    peak = max(calcium) ;
    diastolic = min(calcium) ;
    amplitude = peak - diastolic ;
    value_50 = peak - (0.5 * amplitude) ;
    above_50 = calcium > value_50 ; 
    crossings_Ca = diff(above_50) ;
    up_crossings_Ca = find(crossings_Ca == 1) ;
    down_crossings_Ca = find(crossings_Ca == -1) ;
    if isempty(up_crossings_Ca) || isempty(down_crossings_Ca) %Fallback
        CaD50 = 0 ;
        return ;
    end
    n_pairs_Ca = min(length(up_crossings_Ca), length(down_crossings_Ca)) ;
    CaD50_Diff = time(down_crossings_Ca(1:n_pairs_Ca)) - time(up_crossings_Ca(1:n_pairs_Ca)) ;
    CaD50 = mean(CaD50_Diff) ;
end
