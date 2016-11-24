function main

%  ******* Constants ********
% Set The initial voltage of the membrene to be -60 mili volts
MembraneVoltage = -60; 

% The Time step in our system in mili seconds
TimeStep=0.01; 

% The system time steps array
timeArray=0:TimeStep:100; 

% Calculate the initial values of m,n,h using Aplha & Beta functions
sodiumMGate0 = GetAlphaMGate(MembraneVoltage) / (GetAlphaMGate(MembraneVoltage) + GetBetaMGate(MembraneVoltage));
sodiumHGate0 = GetAlphaHGate(MembraneVoltage) /(GetAlphaHGate(MembraneVoltage) + GetBetaHGate(MembraneVoltage));
potassiumNGate0 = GetAlphaNGate(MembraneVoltage) / (GetAlphaNGate(MembraneVoltage) + GetBetaNGate(MembraneVoltage));

% The values of HH Model in Initial point                
y0=[MembraneVoltage;MembraneVoltage;MembraneVoltage;MembraneVoltage;
    potassiumNGate0;potassiumNGate0;potassiumNGate0;potassiumNGate0;
    sodiumMGate0;sodiumMGate0;sodiumMGate0;sodiumMGate0;
    sodiumHGate0;sodiumHGate0;sodiumHGate0;sodiumHGate0];

tspan = [0,max(timeArray)];


% Call matlab ode45 - runge kutta 4th order method 
% in order to solve Hodgkin-Huxley Model Equation
[timeArray,HHPlot] = ode45(@HodgkinHuxleyEquations,tspan,y0);

figure
subplot(4,1,1);
plot(timeArray, HHPlot(:,1:4));
xlabel('Time (ms)');
ylabel('Voltage (mV)');
legend('node 1', 'node 2', 'node 3', 'node 4');

subplot(4,1,2);
plot(timeArray, HHPlot(:,5:8));
xlabel('Time (ms)');
ylabel('N Varaiable)');
legend('node 1', 'node 2', 'node 3', 'node 4');
subplot(4,1,3);
plot(timeArray, HHPlot(:,9:12));
xlabel('Time (ms)');
ylabel('M Varaiable)');
legend('node 1', 'node 2', 'node 3', 'node 4');
subplot(4,1,4);
xlabel('Time (ms)');
ylabel('H Varaiable)');
plot(timeArray, HHPlot(:,13:16));
legend('node 1', 'node 2', 'node 3', 'node 4');
end

%%Defines the Hodgkin-Huxley Model Equations
function [dydt] = HodgkinHuxleyEquations(t,y)

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

% The resistivity of the axon 
R_a = 35e-4;

% Axon radius
a = 0.5e-4;

% Axon length
dl = 1e-3;

% The resistance of the axon 
R = R_a * dl / (pi * a^2);

% Number of nodes in the axon
x = floor(size(y, 1) / 4);

% Get the external current
I = Current(t);

% Set the params from the function input
V = y(1:x);
n = y(x+ 1:2*x);
m = y(2*x + 1: 3*x);
h = y(3*x + 1: 4*x);

gNa = gbarNa .* m.^3 .* h;
gK = gbarK .* n.^4;
gl = gbarl;
INa = gNa .* (V - ENa);
IK = gK .* (V - EK);
Il = gl .* (V - El);
dydt = zeros(4*x, 1);

%Hodgkin-Huxley Model Equation

%Current equation for all of the parts of the axon
dydt(1:x) = ((-INa - IK - Il)) ./ Cm;

% We add to each 'node' his next 'node' resistance
dydt(1:x-1) = dydt(1:x-1) + dydt(2:x) ./ (Cm*R);

% We add to each 'node' his prev 'node' resistance
dydt(2:x) = dydt(2:x) + dydt(1:x-1) ./ (Cm*R);

% We add the first 'node' the external current
dydt(1) = dydt(1) + (I / Cm);

% The alpha & beta params for each gate
dydt(x+1:2*x) = GetAlphaNGate(V) .* (1 - n) - GetBetaNGate(V) .* n;
dydt(2*x + 1: 3*x) = GetAlphaMGate(V) .* (1 - m) - GetBetaMGate(V) .* m;
dydt(3*x + 1:4*x) = GetAlphaHGate(V) .* (1 - h) - GetBetaHGate(V) .* h;
end

% Alpha value of the m gate 
function alpha = GetAlphaMGate(v) 
 alpha = 0.1 .* (v + 35) ./ (1 - exp(-(v + 35) ./ 10));
end

% Beta value of the m gate 
function beta = GetBetaMGate(v)
beta = 4.0 .* exp(-0.0556 .* (v + 60));
end

% Alpha value of the n gate 
function alpha = GetAlphaNGate(v)
alpha = 0.01 .* (v + 50) ./ (1 -exp(-(v + 50) ./ 10));
end

% Beta value of the n gate 
function beta = GetBetaNGate(v)
beta = 0.125 .* exp(-(v + 60) ./ 80);
end

% Alpha value of the h gate 
function alpha = GetAlphaHGate(v)
alpha = 0.07 .* exp(-0.05 .* (v + 60));
end

% beta value of the h gate 
function beta = GetBetaHGate(v)
beta = 1 ./ (1 + exp(-(0.1) .* (v + 30)));
end