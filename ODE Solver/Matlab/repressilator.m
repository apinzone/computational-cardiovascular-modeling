
%parameter array 
% [alpha, alpha_0, beta, n]

pars = [216, 0, 0.2, 2.0];
%initial conditions 
%  [m1, m2, m3, p1, p2, p3]
y0 = [0 ; 0; 0; 2; 0; 0];

%run for 100s 
t_span = [0 ; 100];

%call solver 
[t, y] = ode15s(@derivatives, t_span, y0, [], pars);

%create plot 
plot(t, y(:,1));
hold on
plot(t, y(:,2));
hold on
plot(t, y(:,3));
hold on 
plot(t, y(:,4));
hold on
plot(t, y(:,5));
hold on
plot(t, y(:,6));
xlabel('Time')
ylabel('mRNA Concentration')

%create derivatives function
function output = derivatives(t, y, p)

    dydt = zeros(size(y));
    %state variables 
    %(m1 , m2 , m3, p1, p2, p3)
    m1 = y(1);
    m2 = y(2);
    m3 = y(3);
    p1 = y(4);
    p2 = y(5);
    p3 = y(6);
    
    %parameters
    alpha = p(1);
    alpha_0 = p(2);
    beta = p(3);
    n = p(4);
    
    %compute derivatives
    dydt(1) = -m1 + alpha/(1+p3.^n) + alpha_0;
    dydt(2) = -m2 + alpha/(1+p1.^n) + alpha_0;
    dydt(3) = -m3 + alpha/(1+p2.^n) + alpha_0;
    dydt(4) = -beta * (p1 - m1);
    dydt(5) = -beta * (p2 - m2);
    dydt(6) = -beta * (p3 - m3);
    
    output = dydt;

end



    