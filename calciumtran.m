%Calcium Transient Model 

%Parameters and Initial Conditions
k_pump = 0.5;
stimulus = 1.0; 
tau = 5.0;
Ca = 0.0;
release = 0.0;
t = 0;
dt = 0.01;
t_end = 50;

times = [t];
calcium_array =[Ca];
release_array = [release];

while t<t_end
    if t < 2.0
        stimulus = 1.0;
    else
        stimulus = 0.0;
    end
    %compute derivatives
    dCa = dCadt(Ca, release, k_pump) * dt;
    drelease = dreleasedt(stimulus, release, tau) * dt;
    %update states
    Ca = Ca + dCa;
    release = release + drelease;
    t = t+dt;
    %Store in arrays
    calcium_array = [calcium_array; Ca];
    release_array = [release_array; release];
    times = [times; t];
end 

plot(times, calcium_array)
hold on
plot(times, release_array)

%Define functions
function dCa = dCadt(Ca, release, k_pump)
    dCa = release - (k_pump * Ca);
end 

function drelease = dreleasedt(stimulus, release, tau)
    drelease = (stimulus - release)/tau;
end 