%Define initial conditions
y0 = [2 , 0] ; %X, Y
t_span = [0, 50] ;

%set up param range
mu_range = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0] ;

%Call solver with parameter sweep
for i = 1:length(mu_range) 
    pars = mu_range(i) ;
    [t, y] = ode15s(@derivatives, t_span, y0, [], pars) ;
    figure
    title(['mu = ' num2str(pars)])
    plot(y(:,1),y(:,2)) ;
    xlabel('x')
    ylabel('y') 

    figure 
    title(['mu = ' num2str(pars)]) ;
    plot(t, y(:,1)) ; %Plot x over time for each mu
    xlabel('Time')
    ylabel('x')
end

function  output = derivatives(t, y, p) 
    derivs = zeros(size(y));
    %Unpack params 
    mu = p(1) ; 
    %State variables
    x = y(1) ;
    y_var = y(2) ; 
    %compute derivatives
    derivs(1) = y_var ;
    derivs(2) = mu * (1 - x.^2) * y_var - x ;
    output = derivs ;
end 