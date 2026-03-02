a = 20 ;
b = 2 ;
c = 5 ; 
%parameters

dt = 0.05;
t_end = 2;
x = c ;
t = 0;
iterations = round(t_end/dt) ;
xall = zeros(iterations, 1) ;

for i = 1:iterations 
	xall(i) = x ;
	dxdt = a - b * x ; %compute derivative
	x = x + dxdt * dt ; %update state variable and multiple derivative by time step at same time
end 

time = dt*(0:iterations-1)' ; %multiply timestep by number of computations to get time for plot
figure
plot(time, xall) %plot time vs. state variables

