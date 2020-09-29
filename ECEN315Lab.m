clf;
clear;
close all;
clc;

Vout=[1.0;2.0;3.0;4.0;5.0;6.0;7.0;8.0];
Iout=[0.13;0.20;0.315;0.47;0.63;0.789;0.995;1.242];
Ra=Vout./Iout;
Ra_Avg=mean(Ra);

V=0.952;
Rs = 220; 
V_t=0.63*V;
te = [0.0067;0.0076];
te = mean(te);
La = te*(Rs+Ra_Avg);

Vb=[1;1.5;2;3;4;5;6]; % voltage 
Ib=[0.0745;0.09;0.115 ;0.173;0.25;0.351;0.45];%Current
T=[0.092,0.0545,0.0358,0.0235,0.0177,0.0171,0.0123];
W=2*pi./T; %angular frequency
% this measure value W, represents the angular velocity, is similar to the
% results in the plotted graph with small differences due to unaccounted
% uncertanties.

Kb = (Vb-Ib*Ra_Avg)./W'
Kb_avg = mean(Kb);
Kt=Kb_avg;
Dm=Kb.*Ib./W';

%change in voltage which .37 lower than the max, from the 0.67 which
%represents the time constant and we are measureing backwards. 
Tm=2*0.260;%s (multiply by 2 because 2 propellers)

Jm = Tm*mean(Dm);
Dm_avg = mean(Dm);

a= Kt/(Jm*La);
b=(Jm*Ra_Avg + Dm_avg*La)/(Jm*La);
c=(Ra_Avg*Dm_avg + Kt^2)/(Jm*La);
sys = tf([a],[1 b c])

step(sys)

hold on
step(2*sys)
step(3*sys)
step(4*sys)
step(5*sys)
step(6*sys)

refline(0,68.3)
refline(0,175.5)
refline(0,267.4)
refline(0,355)
refline(0,432)
refline(0,510)
legend('1V','2V','3V','4V','5V','6V');
title('Step-response of Propeller for Varying Voltages');
xlabel('time t(s)');
ylabel('Rotor angular velocity');


Timeconstant=0.445 % taken off the graphs, where it is pretty constant,
                    % compared to the 0.52 experimental
 
                    
%% lab 3

%Measuring the force using dfferent propeller and scales to measure force
%in grams
%Measuring from 1-6 volts
r=0.13; % m distrance from propeller to axis
Weight_grams = [2,5,9,15,22,31,39,49]; %grams
Force = Weight_grams *0.001*9.81; %newtons
Thrust = Force*r;
Period =[34.6,22.6,17.5,14.1,11.9,10.3,9,8.2];%period in ms
Freq = 1./(Period.*10^-3);
AngularV= 2*pi*(Freq); 

figure
plot(AngularV,Force);
title('Force from Propeller for Varying Angular Velocities');
xlabel('Angular Velocity rad/s');
ylabel('Force N');
%taking the gradient from this we will end up with the value for Kp which
%goes into our final transfer function

%results from excel

A = 0.6904;
B = 0.273;
T_oscl = 0.757;%Period of oscilation
Freq_oscl = 1/T_oscl;
AngFreqOscl = 2*pi*Freq_oscl;
grav = 9.81;%gravitational constant
mass = 0.129; %grams
d= 0.12; % distance to centre of mass of arm from axis

Jp = (mass*d*grav)/(AngFreqOscl^2+B^2); %moment of inertia of the pendulum arm
C =B*2*Jp; % linear damping coefficient

%%

a= Kt/(Jm*La)
b=(Jm*Ra_Avg + Dm_avg*La)/(Jm*La)
c=(Ra_Avg*Dm_avg + Kt^2)/(Jm*La)
Kp = 0.0008; %taken off the graph of excel
d=Jp;
e = C;
f = mass *grav *0.12;
%output = a/(s^2+b*s+c);

figure
sys = tf([a*Kp*r],[d (e+b*d) (f+b*e+c*d) (b*f+e*c) (c*f)])
%sys = tf([a*Kp*r],[1 b c])

step(1*sys)
hold on;
step(2*sys)
step(3*sys)
step(4*sys)
step(5*sys)
step(6*sys)



legend('1V','2V','3V','4V','5V','6V');
title('Angular Displacement Step Response for Varying Voltages ');
xlabel('time t(s)');
ylabel('Angular Displacement (Radians)');

figure
iopzplot(sys)
title('Pole-Zero Map of Transfer Function');





%%
StepV= [1,2,3,4,5,6]; %step voltage in 
StdyAng = [0,2,5,8,10,14]; %steady state angle
GraphAng=[0.0551,0.11,0.166,0.221,0.277,0.331];%radians 
OutputAng =GraphAng.*180/pi;
DifferenceAng =abs(StdyAng-OutputAng)

%%
figure
Finalsys = sys;
ang = 10;
Kp = 1000;
Ki = 1;
control = tf([Kp Ki],[1 0])

M = feedback(Finalsys*0.269*control,1); %multiplying to go from input angle to voltage
%M = feedback(Finalsys*0.269*10,1);
step(ang*M*180/pi)
figure
iopzplot(ang*M*180/pi)
hold on
iopzplot(sys)


