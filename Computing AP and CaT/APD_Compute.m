%Compute APD and CaT from Voltage and Calcium Traces, and Extract Final
%Beat

data = readmatrix('Kernik_values') ;
V = data(:,2) ;
t = data(:, 1) ;
Cai = data(:, 3) ;
BCL = 800 ; %Cycle Length 800ms 
t_end = 60000 ; %Simulation runs for 60 s


%APD Metrics
[APD_30, num_APs_norm_30] = compute_APD(V, t, 30) ;
[APD_60, num_APs_norm_60] = compute_APD(V, t, 60) ;
[APD_90, num_APs_norm_90] = compute_APD(V, t, 90) ;
fprintf('APD at 30 is %.2f\n', APD_30) ;
fprintf('APD at 60 is %.2f\n', APD_60) ;
fprintf('APD at 90 is %.2f\n', APD_90) ;
fprintf('Action Potentials (30) Captured: %.2f\n', num_APs_norm_30) ;
fprintf('Action Potentials (60) Captured: %.2f\n', num_APs_norm_60) ;
fprintf('Action Potentials (90) Captured: %.2f\n', num_APs_norm_90) ;
%CaT Metrics 
[Ca_Peak, Ca_Diastolic, Ca_Amplitude, CaD50] = compute_CaT(Cai, t) ; 
fprintf('Peak Calcium is %.2f\n', Ca_Peak) ;
fprintf('Diastolic Calcium is %.2f\n', Ca_Diastolic) ;
fprintf('Calcium Ampltiude is %.2f\n', Ca_Amplitude) ;
fprintf('CaD50 is %.2f\n', CaD50) ;

%Last Beat Metrics 
last_beat_time = t_end - BCL ; 
last_beat_idx = t >= last_beat_time ; 
V_last_beat = V(last_beat_idx) ; 
t_last_beat = t(last_beat_idx) ; 
Ca_last_beat = Cai(last_beat_idx) ; 
RMP_LB = V_last_beat(1) ; 

%APD Metrics Last Beat 
[APD_30_LB, num_APs_LB_30] = compute_APD(V_last_beat, t_last_beat, 30) ;
[APD_60_LB, num_APs_LB_60] = compute_APD(V_last_beat, t_last_beat, 60) ;
[APD_90_LB, num_APs_LB_90] = compute_APD(V_last_beat, t_last_beat, 90) ;
fprintf('APD at 30 (LB) is %.2f\n', APD_30_LB) ;
fprintf('APD at 60 (LB) is %.2f\n', APD_60_LB) ;
fprintf('APD at 90 (LB) is %.2f\n', APD_90_LB) ;
fprintf('Action Potentials Captured: %.2f\n', num_APs_LB_30) ;
fprintf('RMP (LB) is %.2f\n', RMP_LB) ;

%CaT Metrics Last Beat
[Ca_Peak_LB, Ca_Diastolic_LB, Ca_Amplitude_LB, CaD50_LB] = compute_CaT(Ca_last_beat, t_last_beat) ; 
fprintf('Peak Calcium (LB) is %.2f\n', Ca_Peak_LB) ;
fprintf('Diastolic Calcium (LB) is %.2f\n', Ca_Diastolic_LB) ;
fprintf('Calcium Ampltiude (LB) is %.2f\n', Ca_Amplitude_LB) ;
fprintf('CaD50 (LB) is %.2f\n', CaD50_LB) ;

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
