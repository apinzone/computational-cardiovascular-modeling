%Constants 
%[beta, gamma, delta, l_0, rho, sigma, I_ext] ;
pars = [1, 1, 0.2, 4, 4, 1, 2.5] ; 

%Set initial conditions 
y0 = [2, 1] ; %[ l , lacy]
y0_high = [8, 3] ;
t_span = [0, 20] ; %run for 20 seconds 

%For plotting nullclines must solve derivative equations in terms of one
%variable
%In this case it is easier to put LacY on the Y axis, but also solve for
%both nullclines in terms of LacY

beta = 1;
gamma = 1 ; 
delta = 0.2 ;
l_0 = 4; 
rho = 4 ;
sigma = 1 ; 
l_ext = 2.5 ; 
%create an equation for each nullcline by solving for l and LacY set equal
%to 0
LacY_vals = 0.1:1:10 ;
l_vals = 0.1:1:25 ; 
LacY_nullcline = (delta + rho * (l_vals.^4)./(l_vals.^4 + l_0.^4))/sigma ;
l_nullcline = (beta * l_ext * LacY_vals)/gamma ;
figure
hold on
plot(l_nullcline, LacY_vals, 'b-', 'LineWidth', 2)
plot(l_vals, LacY_nullcline, 'r-', 'LineWidth', 2)
xlabel('l')
ylabel('LacY')
legend('l-nullcline', 'LacY-nullcline')


%solver call
[tv1, yv1] = ode15s(@derivatives, t_span, y0_high, [], pars) ; 
[tv2, yv2] = ode15s(@derivatives, t_span, [3, 1.3], [], pars) ; 
[tv3, yv3] = ode15s(@derivatives, t_span, [3, 1.2], [], pars) ; 
[tv4, yv4] = ode15s(@derivatives, t_span, [2, 1], [], pars) ; 

figure 
hold on
plot(tv1, yv1(:,1), '-r', LineStyle = '--', LineWidth = 2)
plot(tv1, yv1(:,2), '-b', LineStyle = '--', LineWidth = 2)
plot(tv2, yv2(:,1), 'black', LineStyle = '-', LineWidth = 2)
plot(tv2, yv2(:,2), 'blue', LineStyle = '-', LineWidth = 2)
plot(tv3, yv3(:,1), 'magenta', LineStyle = '-.', LineWidth = 2)
plot(tv3, yv3(:,2), 'cyan', LineStyle = '-.', LineWidth = 2)
plot(tv4, yv4(:,1), 'yellow', LineStyle = ':', LineWidth = 2)
plot(tv4, yv4(:,2), 'green', LineStyle = ':', LineWidth = 2)
legend('Lactose [8,3]','LacY [8,3]', 'Lactose [3,1.3]','LacY [3,1.3]','Lactose [3, 1.2]','LacY [3, 1.2]', 'Lactose [2,1]','LacY [2,1]') 

%param sweep for l_ext 
l_ext_range = 1:0.5:7 ; 
LacY_high_final = [] ;
LacY_base_final = [] ;
for i = 1:length(l_ext_range)
    pars = [1, 1, 0.2, 4, 4, 1, l_ext_range(i)] ;
    %l = 8, LacY = 3 
    [t1,y1] = ode15s(@derivatives, t_span, y0_high, [], pars) ;
    %l = 2, LacY = 1 
    [t2,y2] = ode15s(@derivatives, t_span, y0, [], pars) ;
    %Get final values of LacY over l_ext range
    LacY_high_final(i) = y1(end, 2) ;
    LacY_base_final(i) = y2(end, 2) ;
end

figure 
hold on 
plot(l_ext_range,LacY_base_final, '-r', 'LineWidth', 2)
plot(l_ext_range,LacY_high_final, '-b', 'LineWidth', 2)
legend('High Initial Conditions [8, 3]', 'Base  Initial Conditions [2, 1]')
xlabel('l_ext')
ylabel('LacY')

function output = derivatives(t, y, p)
    derivs = zeros(size(y)) ;
    %unpack params
    beta = p(1) ;
    gamma = p(2) ; 
    delta = p(3) ; 
    l_0 = p(4) ; 
    p_var = p(5) ;
    sigma = p(6) ; 
    l_ext = p(7) ; 
    %get state variables 
    l = y(1) ; 
    LacY = y(2) ; 
    %compute derivatives 
    derivs(1) = beta * l_ext * LacY - gamma * l ; %lactose
    derivs(2) = delta * p_var * (l.^4/(l.^4 + l_0.^4)) - sigma * LacY ; %LacY
    output = derivs; %return derivatives 
end
