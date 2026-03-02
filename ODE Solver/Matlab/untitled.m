%Gene Bistability Toggle Switch Model, Gardner et al. 2000

%parmeter sweep setup 
n_range = linspace(1, 5, 20) ; %at y0 = [1.5, 1.5]
alpha_range = linspace(1, 4, 30); 
base_pars = [2.5, 1, 2] ; %alpha, beta, n 

%Initial conditions and timespan
y0_B_high = [0.5, 3.0] ;
y0_A_high = [3.0, 0.5] ;
y0 = [1.5, 1.5] ;
t_span = [0 , 50] ; 

%Solver Call Basic (Three Different y0) 
[t1, y1] = ode15s(@derivatives, t_span, y0_B_high , [], base_pars) ;
[t2, y2] = ode15s(@derivatives, t_span, y0_A_high , [], base_pars) ;
[t3, y3] = ode15s(@derivatives, t_span, y0 , [], base_pars) ;

figure 
hold on
%Plot with lower opacity 
plot(t1, y1(:,1), 'Color', [0, 0, 1, 0.6], 'LineWidth', 2, 'LineStyle', '-');  % A at High B 
plot(t1, y1(:,2), 'Color', [1, 0, 0, 0.6], 'LineWidth', 2, 'LineStyle', '--');  % B at High B 
plot(t2, y2(:,1), 'Color', [0, 1, 0, 0.6], 'LineWidth', 2, 'LineStyle', '-');  % A at High A 
plot(t2, y2(:,2), 'Color', [1, 1, 0, 0.6], 'LineWidth', 2, 'LineStyle', '--');  % B at High A 
plot(t3, y3(:,1), 'Color', [0, 1, 1, 0.6], 'LineWidth', 2, 'LineStyle', '-');  % A at moderate 
plot(t3, y3(:,2), 'Color', [1, 0, 1, 0.6], 'LineWidth', 2, 'LineStyle', '--');  % B at moderate 
legend('A at High B', 'B at High B', 'A at High A', ...
    'B at High A', 'A at moderate', 'B at moderate') ;
xlabel('Time')
ylabel('Protein Concentration')

%Solver Call with parameter sweep for hill coefficient n 
A_final_array = zeros(size(n_range)) ;
B_final_array = zeros(size(n_range)) ;
for i = 1:length(n_range)
    pars = [2.5, 1, n_range(i)] ; 
    [t4, y4] = ode15s(@derivatives, t_span, y0, [], pars) ;
    A_column = y4(:, 1) ;
    A_final_array(i) = A_column(end) ;
    B_column = y4(:, 2) ;
    B_final_array(i) = B_column(end) ;
end

%Plot hill parameter sweep
figure 
hold on
plot(n_range, A_final_array, '-r', 'LineStyle', '-')
plot(n_range, B_final_array, '-b', 'LineStyle', '--')
title('A and B concentrations across hill coefficients (1 to 5)')
legend('A', 'B')
xlabel('n')
ylabel('Protein Concentration')

%Solver Call with parameter sweep for alpha
bifurcation_A_array_B_high = zeros(size(alpha_range)) ;
bifurcation_B_array_B_high = zeros(size(alpha_range)) ;
bifurcation_A_array_A_high = zeros(size(alpha_range)) ;
bifurcation_B_array_A_high = zeros(size(alpha_range)) ;

for i = 1:length(alpha_range)
    pars = [alpha_range(i), 1, 2] ; 
    %B high sim for param sweep alpha
    [t5, y5] = ode15s(@derivatives, t_span, y0_B_high, [], pars) ;
    %A high sim for param sweep alpha
    [t6, y6] = ode15s(@derivatives, t_span, y0_A_high, [], pars) ;
    %Get final values for B high simulations
    A_column_B_high = y5(:, 1) ;
    bifurcation_A_array_B_high(i) = A_column_B_high(end) ;
    B_column_B_high = y5(:, 2) ;
    bifurcation_B_array_B_high(i) = B_column_B_high(end) ;
    %Get final values for A high simulations
    A_column_A_high = y6(:, 1) ;
    bifurcation_A_array_A_high(i) = A_column_A_high(end) ;
    B_column_A_high = y6(:, 2) ;
    bifurcation_B_array_A_high(i) = B_column_A_high(end) ;
end

%Plot alpha bifurcation 
figure 
hold on
%Plot
plot(alpha_range, bifurcation_A_array_B_high, 'Color', [1, 0, 0, 0.7], 'LineWidth', 2.5, 'LineStyle', '-');   % A at B-high (red solid)
plot(alpha_range, bifurcation_A_array_A_high, 'Color', [0, 0, 1, 0.7], 'LineWidth', 2.5, 'LineStyle', '--');  % A at A-high (blue dashed)
plot(alpha_range, bifurcation_B_array_B_high, 'Color', [0, 0.8, 0, 0.7], 'LineWidth', 2.5, 'LineStyle', '-'); % B at B-high (green solid)
plot(alpha_range, bifurcation_B_array_A_high, 'Color', [0, 0, 0, 0.7], 'LineWidth', 2.5, 'LineStyle', '--');  % B at A-high (black dashed)
legend('A with high initial B', 'A with high initial A', ...
    'B with high initial B', 'B with high initial A')
xlabel('Alpha Value')
ylabel('Protein Concentration') 

%create derivatives function
function output = derivatives(t, y, p)
    %blank derivatives array 
    derivs = zeros(size(y)) ;
    %Unpack params
    alpha = p(1) ; 
    beta = p(2) ; 
    n = p(3) ; 
    %Get state variables
    A = y(1) ;
    B = y(2) ;
    %compute derivatives and updated derivs
    derivs(1) = alpha / (1 + B.^n) - A ;%dA/Dt
    derivs(2) = alpha / (1 + A.^n) - B ; %dBdt
    %return array of derivatives 
    output = derivs ; 
end 
