% Main File to generate Baseline Model Figures 10-11
close all; clear; clc

load ICs_baseline
load baseline_parameter_inputs
baseline_parameter_inputs(2)  = 0.7 ;   % IKr
baseline_parameter_inputs(5)  = 1.2 ;   % ICaL
baseline_parameter_inputs(8)  = 1.3 ;   % If
baseline_parameter_inputs(9)  = 0.5 ;   % NCX
baseline_parameter_inputs(10) = 0.5 ;   % SERCA
baseline_parameter_inputs(11) = 1.5 ;   % RyR
baseline_parameter_inputs(15) = 1.3 ;   % IbCa
baseline_parameter_inputs(16) = 0.6 ;   % PMCA

%% Run iPSC_function
options = odeset('MaxStep',1,'InitialStep',2e-2);
run_time=10.6e4;
[Time, values] = ode15s(@ipsc_function,[0, run_time],Y_init, options, baseline_parameter_inputs);
[Time2, values2] = ode15s(@ipsc_function_2000,[0, run_time],Y_init, options, baseline_parameter_inputs);
last_6_idx = Time >= 10e4 ; 
last_6_idx_t2 = Time2 >= 10e4 ;
time_trimmed = Time(last_6_idx) ; 
time_trimmed_2 = Time2(last_6_idx_t2) ; 
Cai=values(last_6_idx,3);
Vm=values(last_6_idx,1);
Cai2=values2(last_6_idx_t2,3);
Vm2=values2(last_6_idx_t2,1);

%Write Calcium and Voltage Time Series to File 
values_table_1000 = table(time_trimmed, Vm, Cai) ;
writetable(values_table_1000, 'Kernik_values_1000.csv')
values_table_2000 = table(time_trimmed_2, Vm2, Cai2) ;
writetable(values_table_2000, 'Kernik_values_2000.csv')

%% Calculate select current traces:
INaCa = zeros(size(Time));
IpCa = zeros(size(Time));
Iup = zeros(size(Time));
for i= 1:size(values,1)
    [~, update_step_i] =  ipsc_function(Time(i), values(i,:),  baseline_parameter_inputs);    
    INaCa(i) = update_step_i(8);
    IpCa(i) = update_step_i(9);
    Iup(i) = update_step_i(14);
end

%% Figure 10A & 10C: Calcium Flux analysis and Calcium Transient Trace
ca_analysis( Time, Iup, INaCa, IpCa, Cai )

%% Figure 11A: action potential trace for baseline model 
figure,set(gcf,'color','w')
plot(Time, Vm,'Color', [.8 0 .18]);
set(gca,'box','off','tickdir','out')
ylabel('Voltage (mV)');
xlabel('Time (ms)')


