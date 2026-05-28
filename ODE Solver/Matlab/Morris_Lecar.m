%Morris Lecar 

%Fixed Params 
%           [C, g_Ca, g_K, g_L, V_Ca, V_K, V_L, V1, V2, V3, V4, phi]
base_pars = [20, 4.4 , 8, 2, 120, -84, -60, -1.2, 18, 2, 30, 0.04] ; 

%Initial conditions
y0 = [-60, 0] ; % V0, w0
t_span = [0, 500] ; %ms 

I_app_range = [20, 50, 80] ;
Voltage_max = zeros(size(I_app_range)) ;
Voltage_min = zeros(size(I_app_range)) ;
%Solver call for range of I_apps 
for i = 1:length(I_app_range)
    I_app = I_app_range(i) ;
    pars = [base_pars, I_app] ;
    [t, y] = ode15s(@derivatives, t_span, y0, [], pars) ;
    
    %Get idx for last 50%
    settled_idx = round(length(t) / 2); 
    V_settled = y(settled_idx:end, 1) ;
    Voltage_max(i) = max(V_settled) ;
    Voltage_min(i) = min(V_settled) ;

    %Plot
    figure 
    plot(t, y(:,1)) ;

    figure
    plot(y(:,1), y(:,2)) ;

end


%Bifurcation plot 
figure 
hold on
plot(I_app_range, Voltage_max, '-r') ;
plot(I_app_range, Voltage_min, '-b') ;
legend('Max Voltage', 'Min Voltage') 
xlabel('I_app')
ylabel('Voltage (mV)') 

function output = derivatives(t, y, p) 
    %Create derivatives array
    derivs = zeros(size(y)) ;
    %unpack pars 
    C = p(1) ;
    g_Ca = p(2) ;
    g_K = p(3) ; 
    g_L = p(4) ; 
    V_Ca = p(5) ; 
    V_K = p(6) ; 
    V_L = p(7) ;
    V1 = p(8) ;
    V2 = p(9) ; 
    V3 = p(10) ;
    V4 = p(11) ;
    phi = p(12) ;
    I_app = p(13) ; 

    %Get state variables
    V = y(1) ;
    w = y(2) ; 
    
    %Relevant calculations
    m_inf = 0.5 * (1 + tanh((V - V1) / V2)) ;
    w_inf = 0.5 * (1 + tanh((V - V3)/ V4)) ;
    t_w = 1 / cosh((V - V3) / (2 * V4)) ;

    %Compute derivatives
    derivs(1) = (I_app - g_Ca * m_inf * (V - V_Ca) - g_K * w * (V - V_K) - g_L * (V - V_L)) / C ; %dVdt 
    derivs(2) = phi * (w_inf - w) / t_w ;

    %Return Output 
    output = derivs ;
end 

   