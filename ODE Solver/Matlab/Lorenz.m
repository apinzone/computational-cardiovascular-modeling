
%params 
pars = [10, 28, 8/3]; %sigma, rho, beta 

%initial conditions
t_span = [0, 50] ; %run for 50s 
y0 = [1 , 1, 1] ; %[x , y, z]

%call solver 
[t, y] = ode15s(@derivatives, t_span, y0, [], pars) ;
%plot 
figure
hold on
plot(t, y(:,1), 'r') ;
plot(t, y(:,2), 'b') ;
plot(t, y(:,3), 'g') ;
legend('x', 'y', 'z')

figure
plot3(y(:,1), y(:,2), y(:,3)) ;

function output = derivatives(t, y, p)
    derivs = zeros(size(y)) ;
    %unpack params
    sigma = p(1) ;
    rho = p(2) ;
    beta = p(3) ; 
    %state variables
    x = y(1) ;
    y_var = y(2) ;
    z = y(3) ;
    %compute derivatives
    derivs(1) = sigma * (y_var - x) ;
    derivs(2) = x * (rho - z) - y_var ;
    derivs(3) = x*y_var - beta * z ;
    %return derivatives 
    output = derivs ;
end

