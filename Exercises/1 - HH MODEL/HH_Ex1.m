%Rebuild Hodgkin-Huxley Model
%Initial Conditions and parameters
t_span = [0, 25] ; %Run for 25 ms
y0 = [-60, 0.3177, 0.0529, 0.5961] ; %V, n, m, h
I_stim_range = linspace(0.05, 1.8, 15) ; %Range of 15 values from 0.05 to 1.8 to vary stimulus current
pulse_start = 2 ; %Start pulse at 2 ms 
pulse_end_range = [3, 3.25, 3.5, 3.75, 4, 5, 6, 7, 9, 12] ; %1ms, 1.25ms, 1.5ms, 1.75ms, 2ms, 3ms, 4ms, 5ms, 7 ms, 10 ms
figure %Start figure before loop for subplot

%Solver Call and Variable Extraction, test different values for I_stim 
for i = 1:length(I_stim_range)
    I_stim = I_stim_range(i) ;
    pars = [0.01, 55.17, -72.14, -49.42, 1.2, 0.36, 0.003, I_stim, 0, 0] ; 
    %1.C_m, 2.E_Na, 3.E_K, 4.E_I, 5.g_Na, 6.g_K, 7.g_leak, 8.stim current amp 9. pulse_start, 10. pulse_end (0 for 9 and 10 if constant stimulus)
    [t, y] = ode15s(@derivatives, t_span, y0, [], pars) ;
    %Plot Voltage over time for each I_stim value - Figure 1. 
    subplot(3, 5, i)
    plot(t, y(:,1)) ;
    xlabel('Time (ms)')
    ylabel('Voltage (mv)')
    ylim([-80 60])
    xlim([0, 25])
    yticks(-80:20:60)
    xticks(0:5:25)
    title('I stim = ', I_stim)
end

figure
% Solver call for testing different pulse durations
for i = 1:length(pulse_end_range)
    %Use "best" stimulus current amplitude from sweep (0.175)
    pars_pulse = [0.01, 55.17, -72.14, -49.42, 1.2, 0.36, 0.003, 0.175, pulse_start, pulse_end_range(i)] ;
    %1.C_m, 2.E_Na, 3.E_K, 4.E_I, 5.g_Na, 6.g_K, 7.g_leak, 8.stim current amp, 9. pulse start time, 10. pulse end time
    pulse_duration = pars_pulse(10) - pars_pulse(9) ;
    [t_pulse, y_pulse] = ode15s(@derivatives, t_span, y0, [], pars_pulse) ;
    %Plot Voltage over time for each I_stim value - Figure 2.
    subplot(2, 5, i)
    plot(t_pulse, y_pulse(:,1)) ;
    xlabel('Time (ms)')
    ylabel('Voltage (mv)')
    ylim([-80 60])
    xlim([0, 25])
    yticks(-80:20:60)
    xticks(0:5:25)
    title('Pulse Duration = ', pulse_duration)
end

%Best parameters solver call for final figure
pars_best = [0.01, 55.17, -72.14, -49.42, 1.2, 0.36, 0.003, 0.175, 2, 3] ; %'Best' pulse duration of 1.5ms, I_stim of 0.175
    %1.C_m, 2.E_Na, 3.E_K, 4.E_I, 5.g_Na, 6.g_K, 7.g_leak, 8.stim current amp, 9. pulse start time, 10. pulse end time
[t_best, y_best] = ode15s(@derivatives, t_span, y0,[], pars_best) ;
V_best = y_best(:, 1) ;
%Plot Voltage over time - Figure 3.
figure
plot(t_best, V_best) 
xlabel('Time (ms)')
ylabel('Voltage (mv)')
title('Action Potential: I stim = 0.175, Pulse Duration = 1 ms')
ylim([-80 60])
xlim([0, 25])
yticks(-80:20:60)
xticks(0:5:25)

%Derivatives Function
function output = derivatives(t,y,p) 
    %blank deriv array size of state variables
    derivs = zeros(size(y)) ;
    %Unpack state variables 
    V = y(1) ;
    n = y(2) ;
    m = y(3) ;
    h = y(4) ; 
    %Unpack pars 
    C_m = p(1) ; 
    E_Na = p(2) ;
    E_K = p(3) ;
    E_leak = p(4) ;
    g_Na = p(5) ;
    g_K = p(6) ; 
    g_leak = p(7) ;
    t_pulse_start = p(9) ;
    t_pulse_end = p(10);
    if t_pulse_start == 0 && t_pulse_end == 0 %option for constant stimulus application
        I_stim = p(8) ;
    else 
        if t >= t_pulse_start && t <= t_pulse_end %option for pulse 
            I_stim = p(8) ;
        else
            I_stim = 0 ;
        end
    end

    %Calculate rate constants 
    alpha_n = 0.01 * (V + 50) / (1 - exp(-(V + 50)/10)) ;
    beta_n = 0.125 * exp((-(V+60))/80) ;
    alpha_m = 0.1 * (V + 35) / (1 - exp(-(V + 35)/10)) ;
    beta_m = 4 * exp(-0.0556*(V + 60)) ;
    alpha_h = 0.07 * exp(-0.05 * (V + 60)) ;
    beta_h = 1 / (1 + exp(-0.1 * (V + 30))) ;
    %Calculate Currents
    I_Na = g_Na * m.^3 * h * (V - E_Na) ;
    I_K = g_K * n.^4 * (V - E_K) ;
    I_leak = g_leak * (V - E_leak) ;
    %Compute Derivatives
    derivs(1) = (1/C_m) * (I_stim - I_Na - I_K - I_leak) ; %dVdt
    derivs(2) = alpha_n * (1 - n) - beta_n * n ; %dndt 
    derivs(3) = alpha_m * (1 - m) - beta_m * m ; %dmdt
    derivs(4) = alpha_h * (1 - h) - beta_h * h ; %dhdt 
    %Return 
    output = derivs ;
end
