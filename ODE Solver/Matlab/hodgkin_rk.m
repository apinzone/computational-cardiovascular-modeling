
%Initial Conditions
y0 = [-65.0; 0.05; 0.6; 0.32]; %V, m , h, n
t_span = [0; 50]; %0 - 50s. span

%parameters 
Cm = 1.0; %membrane capacitance
gNa = 120; %max sodium conductance
gK = 36; %max potassium conductance
gL = 0.3; %leak conductance
ENa = 50; %sodium reversal
EK = -77.0; %potassium reversal
EL = -54.4; %leak reversal
I_ext = 10; %ext stimulus current


%call solver
params = [Cm, gNa, gK, gL, ENa, EK, EL, I_ext];
[t,y] = ode15s(@derivatives,t_span,y0, [], params);
%plot
time_array = t;
y_array = y;
plot(t, y(:,1));  
%Define Derivative function
function output = derivatives(t, y, p)
    %unpack parameters
    Cm = p(1);
    gNa = p(2);
    gK = p(3);
    gL = p(4);
    ENa = p(5);
    EK = p(6);
    EL = p(7);
    I_ext = p(8);

    %get state variables
    V = y(1);
    m = y(2);
    h = y(3);
    n = y(4);
    %compute currents and gates 
    INa = gNa * m.^3 * h * (V - ENa);
    IK = gK * n.^4 * (V - EK);
    IL = gL * (V -EL);
    alpha_m = 0.1 * (V + 40) / ( 1 - exp(-(V + 40)/10));
    beta_m = 4.0 * exp(-(V + 65) / 18);
    alpha_h = 0.07 * exp(-(V + 65) / 20);
    beta_h = 1.0 / (1 + exp(-(V + 35)/10));
    alpha_n = 0.01 * (V + 55) / (1 - exp(-(V + 55)/10));
    beta_n = 0.125 * exp(-(V + 65)/80);
    %compute derivatives
    dVdt = (1/Cm) * (I_ext - INa - IK - IL);
    dmdt = alpha_m * (1 - m)- beta_m * m;
    dhdt = alpha_h * (1-h) - beta_h * h;
    dndt = alpha_n * (1-n) - beta_n * n;
    %update arrays
    output = [dVdt; dmdt; dhdt; dndt];
end

