function main

%  ******* Constants ********
% Set The initial voltage of the membrene to be -60 mili volts
MembraneVoltage = -60; 

% The Time step in our system in mili seconds
TimeStep=0.01; 

% The system time steps array
timeArray=0:TimeStep:25; 

% Calculate the initial values of m,n,h using Aplha & Beta functions
sodiumMGate0 = GetAlphaMGate(MembraneVoltage) / (GetAlphaMGate(MembraneVoltage) + GetBetaMGate(MembraneVoltage));
sodiumHGate0 = GetAlphaHGate(MembraneVoltage) /(GetAlphaHGate(MembraneVoltage) + GetBetaHGate(MembraneVoltage));
potassiumNGate0 = GetAlphaNGate(MembraneVoltage) / (GetAlphaNGate(MembraneVoltage) + GetBetaNGate(MembraneVoltage));


% The values of HH Model in Initial point                
y0=[MembraneVoltage; potassiumNGate0; sodiumMGate0; sodiumHGate0];

tspan = [0,max(timeArray)];


% Call matlab ode45 - runge kutta 4th order method 
% in order to solve Hodgkin-Huxley Model Equation
[timeArray,HHPlot] = ode113(@HodgkinHuxleyEquations,tspan,y0);

% Split the ode output matrix
VolatgePlot = HHPlot(:,1);
nChangePlot = HHPlot(:,2);
mChangePlot = HHPlot(:,3);
hChangePlot = HHPlot(:,4);

%Plot the functions
figure
title('Hodgkin-Huxley Model');

subplot(3,1,1);
plot(timeArray,VolatgePlot);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('Voltage Change - Hodgkin-Huxley Model');

subplot(3,1,2);
plot(timeArray, nChangePlot,'b' , timeArray, mChangePlot,'g', timeArray, hChangePlot, 'r');
ylabel('Gaining Variables')
xlabel('Time (ms)')
legend('n', 'm', 'h');

subplot(3,1,3);
x = 0:1:25;
y =  0.2*sigmf(x,[10000 1]) -0.2*sigmf(x,[10000 2]) + 0.2*sigmf(x,[10000 10]) - 0.2*sigmf(x,[10000 11]) + 0.2*sigmf(x,[10000 14]) - 0.2*sigmf(x,[10000 15]);
plot(x,y);
xlabel('Time (ms)');
ylabel('Current');
ylim([0.05 0.35]);
xlim([0,25]);
end

%%Defines the Hodgkin-Huxley Model Equations
function dydt = HodgkinHuxleyEquations(t,y)

%  Na reversal potential in mili volts
ENa=55.17; 

% K reversal potential in mili volts
EK=-72.14;

% Leakage reversal potential in mili volts
El=-49.42;

%  Na conductance in mS/cm^2
gbarNa=1.2;

% K conductance in mS/cm^2 
gbarK=0.36;

% Leakage conductance in mS/cm^2
gbarl=0.003;

%The Membrane Capacitance in uF/cm^2
Cm = 0.01; 


%{
%Change the current after 1/2 of the time period
if (t < 8)
I = 0.1;
% The current value
else
I = 0.3; 
end;
% I mu + sigma*randn
%}
I = 0.2*sigmf(t,[10000 1]) -0.2*sigmf(t,[10000 2]) + 0.2*sigmf(t,[10000 10]) - 0.2*sigmf(t,[10000 11]) + 0.2*sigmf(t,[10000 14]) - 0.2*sigmf(t,[10000 15]);

%I =  0.1 + 0.2*sigmf(t,[10000 13]);
%{
%Normal distribiution calc for the current
 mu = 0.03;
 sigma = 0.9;
% 
 I = mu + sigma * randn;
%end
%}

% Set the params from the function input
V = y(1);
n = y(2);
m = y(3);
h = y(4);

gNa = gbarNa * m^3 * h;
gK = gbarK * n^4;
gl = gbarl;
INa = gNa * (V - ENa);
IK = gK * (V - EK);
Il = gl * (V - El);

%Hodgkin-Huxley Model Equation
dydt = [((1 / Cm) * (I - (INa + IK + Il)));
        GetAlphaNGate(V) * (1 - n) - GetBetaNGate(V) * n;
        GetAlphaMGate(V) * (1 - m) - GetBetaMGate(V) * m;
        GetAlphaHGate(V) * (1 - h) - GetBetaHGate(V) * h];
end

% Alpha value of the m gate 
function alpha = GetAlphaMGate(v) 
 alpha = 0.1 * (v + 35) / (1 - exp(-(v + 35) / 10));
end

% Beta value of the m gate 
function beta = GetBetaMGate(v)
beta = 4.0 * exp(-0.0556 * (v + 60));
end

% Alpha value of the n gate 
function alpha = GetAlphaNGate(v)
alpha = 0.01 * (v + 50) / (1 -exp(-(v + 50) / 10));
end

% Beta value of the n gate 
function beta = GetBetaNGate(v)
beta = 0.125 * exp(-(v + 60) / 80);
end

% Alpha value of the h gate 
function alpha = GetAlphaHGate(v)
alpha = 0.07 * exp(-0.05 * (v + 60));
end

% beta value of the h gate 
function beta = GetBetaHGate(v)
beta = 1 / (1 + exp(-(0.1) * (v + 30)));
end

