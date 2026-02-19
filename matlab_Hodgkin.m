%Parameters
Cm = 1.0; %membrane capacitance
gNa = 120; %max sodium conductance
gK = 36; %max potassium conductance
gL = 0.3; %leak conductance
ENa = 50; %sodium reversal
EK = -77.0; %potassium reversal
EL = -54.4; %leak reversal
I_ext = 10; %ext stimulus current

%Initial Conditions
V = -65.0;
m = 0.05;
h = 0.6;
n = 0.32;
t=0;
dt = 0.01;
t_end = 50; 

%Create Arrays 
time_array = t;
voltage_array = V; 
m_array = m;
h_array = h;
n_array = n;

%Solver Loop
while t<t_end   
%Compute Derivatives
    dV = dVdt(V,Cm,m, h, n, gNa, gK, gL, ENa, EK, EL, I_ext) * dt;
    dm = dmdt(V, m) * dt;
    dh = dhdt(V, h) * dt;
    dn = dndt(V, n) * dt;
%Update State Variables
    V = V + dV;
    m = m + dm;
    h = h + dh;
    n = n + dn;
    t = t + dt;
%Store update variables in arrays 
    voltage_array = [voltage_array ; V];
    m_array = [m_array ; m];
    h_array = [h_array ; h];
    n_array = [n_array; n];
    time_array = [time_array; t];
end

figure(1)
plot(time_array, voltage_array);
xlabel('Time (s)');
ylabel('Voltage (mV)');

figure(2)
plot(time_array, m_array);
hold on
plot(time_array, h_array);
hold on
plot(time_array, n_array);
xlabel('Time (s)');
ylabel('Gating Activity');
legend('m', 'h', 'n');

%Define Equations
function dV = dVdt(V,Cm,m, h, n, gNa, gK, gL, ENa, EK, EL, I_ext)
    INa = gNa * m.^3 * h * (V - ENa);
    IK = gK * n.^4 * (V - EK);
    IL = gL * (V -EL);
    dV = (1/Cm) * (I_ext - INa - IK - IL);
end 

function dm = dmdt(V, m)
    alpha_m = 0.1 * (V + 40) / ( 1 - exp(-(V + 40)/10));
    beta_m = 4.0 * exp(-(V + 65) / 18);
    dm = alpha_m * (1 - m)- beta_m * m;
end

function dh = dhdt(V, h)
    alpha_h = 0.07 * exp(-(V + 65) / 20);
    beta_h = 1.0 / (1 + exp(-(V + 35)/10));
    dh = alpha_h * (1-h) - beta_h * h;
end

function dn = dndt(V, n)
    alpha_n = 0.01 * (V + 55) / (1 - exp(-(V + 55)/10));
    beta_n = 0.125 * exp(-(V + 65)/80);
    dn = alpha_n * (1-n) - beta_n * n;
end