%ODE course my way 
a = 20 ;
b = 2 ;
c = 5 ; 
%parameters

dt = 0.01;

%initial conditions
t = 0;
t_end = 2;
x = c ;

x_array = [x]; 
times = [t] ;

while t<t_end
   %compute derivatives 
   dx = (a - b * x ) * dt ;
   %update state variables
   x = x + dx; 
   t = t + dt ; 
   %store in array 
   x_array = [x_array ; x] ;
   times = [times ; t] ;
end 

plot(times, x_array)