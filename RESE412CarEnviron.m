
clf;
clear;
close all;
clc;

%% Processing Data
load('sodaprodata.mat')
SodaProRaw = sodaprodata;
OctoberWindspeed = zeros(24,30);
sodaprodata.Hour=grp2idx(sodaprodata.Hour);
sodaprodata = cell2mat(table2cell(sodaprodata));

%Setting Variables
N = length(sodaprodata); %Total amount of rows
Month = 10;
Info = 0;% what collumn to take data from
Increment = 1;
SolarData = [];
T = 1:24;
%only observing octber
for a = 1:N
   if sodaprodata(a,2) == Month
      SolarData(Increment,:) = sodaprodata(a,:);
      Increment = Increment +1;
    
   end
end
OctoberRaw = SolarData;

Temp = [];
inc = 1;
for i = 1:744:length(SolarData)
    Temp(:,:,inc) = SolarData(i:i+743,:);
    inc = inc+1;
end
SolarData = Temp;




%% Organise the data
 
DataAvg = mean(Temp,3);
DataVar = var(Temp,0,3);

Temp1 = (1:12)';
Temp2 = (13:24)';
Temporary = [];
for i = 1:24:length(DataAvg)
    DataAvg(i:i+11,4) = Temp2;
    DataVar(i:i+11,4) = Temp2;
    DataAvg(i+12:i+23,4) = Temp1;
    DataVar(i+12:i+23,4) = Temp1;   
end

for i = 1:24:length(OctoberRaw)
    OctoberRaw(i:i+11,4) = Temp2;
    OctoberRaw(i+12:i+23,4) = Temp1; 
end



%% Plotting WindSpeeds
OctWindSpd = zeros(24,1);
count = 0;
for a = 1:length(DataAvg)   
        for b = 1:24 %iterate each hour
            if DataAvg(a,4) == b
                OctWindSpd(b,1) =  OctWindSpd(b,1) + DataAvg(a,8)/30;%Windspeed 
                count = count +1;
            end           
        end   
end


figure
plot(T,OctWindSpd)
ylabel('WindSpeed m/s') 
ylim([5 8])
grid on
xlim([0 23])
xticks([0 2 4 6 8 10 12 14 16 18 20 22 23])
xlabel('Time Hours') 
title(['Average Windspeed over a Day in October'])


%% Plotting Solar Irradiance
OctIrradiance = zeros(24,1);
OctIrradianceVar = zeros(24,1);
count = 0;
for a = 1:length(DataAvg)   
        for b = 1:24 %iterate each hour
            if DataAvg(a,4) == b
                OctIrradianceVar(b,1) =  OctIrradianceVar(b,1) + DataVar(a,13)/30;
                OctIrradiance(b,1) =  OctIrradiance(b,1) + DataAvg(a,13)/30;%Windspeed 
                count = count +1;
            end           
        end   
end


%figure
%errorbar(T,OctIrradiance,OctIrradianceVar)
%plot(T,OctIrradiance)

figure
hold on
plot(T,OctIrradiance)
scatter(DataAvg(:,4),DataAvg(:,13))
grid on
ylabel('Irradiance W/m^2') 
%ylim([0 18])
%xlim([0 23])
xticks([0 2 4 6 8 10 12 14 16 18 20 22 23])
xlabel('Time Hours') 
title(['Average Irradiance over a Day in October'])
legend('Average Irradiance','Average Points of 2000-19')


figure
plot(T,OctIrradianceVar)
grid on
ylabel('Irradiance Variance') 
%ylim([0 18])
%xlim([0 23])
xticks([0 2 4 6 8 10 12 14 16 18 20 22 23])
xlabel('Time Hours') 
title(['Irradiance Variance over a Day in October for last 19 years'])
%legend('Average Irradiance')

figure
hold on
scatter(OctoberRaw(:,4),OctoberRaw(:,13),'x')
plot(T,OctIrradiance,'LineWidth',2.0)
xlabel('Time Hours')
ylabel('Solar Irradiance W/m^2')
title(['Solar Irradiance Data Points over 2000-2019 with Averagee'])
legend('Data Points','Average Irradiance')



%% determining the total power available
%Polycrystaline efficiency = 14-16%
F = griddedInterpolant(T,OctIrradiance);
fun = @(t) F(t);
PotentialPower_kWh = integral(fun,T(10), T(16))/1000
Efficiency = 15
EffectivePower_kWh = PotentialPower_kWh * Efficiency/100



%% Plotting Temperature Averages
OctTemp = zeros(24,1);
OctTempVar = zeros(24,1);
count = 0;
for a = 1:length(DataAvg)   
        for b = 1:24 %iterate each hour
            if DataAvg(a,4) == b
                OctTempVar(b,1) =  OctTempVar(b,1) + DataVar(a,5)/31;
                OctTemp(b,1) =  OctTemp(b,1) + DataAvg(a,5)/31;%Windspeed 
                count = count +1;
            end           
        end   
end



OctTemp = OctTemp-273;%change from kelvin to celsius
figure
hold on
errorbar(T,OctTemp,OctTempVar)
plot(T,OctTemp,'r')
grid on
ylabel('Temperature C') 
ylim([0 20])
xlim([0 23])
xticks([0 2 4 6 8 10 12 14 16 18 20 22 23])
xlabel('Time Hours') 
title(['Average Temperature over a Day in October'])
legend('Temperature Error bars of Averages','Average Temperature')

%%
Windspeed = 0:20;
WindPower = [0,0,0,2.5,5,10,17,25,34,41,50,56,60.5,64,69,71.5,73,74,75,75,75];
figure
hold on
plot(Windspeed,WindPower)
grid on
xlabel('Wind Speed m/s') 
ylabel('Power W') 
title(['Windspeed and Turbine power relationship'])
W = OctWindSpd;
%Power = -3*10^-5*6* + 0.0019x5 - 0.0447x4 + 0.412x3 - 0.8385x2 + 0.4494x + 0.022;


%% Data Analysis of RC Car Speed
%% Proper data test 1

SampleSpeed = 0.02025;
Distance = 1.45;
 load('test4.mat')
test4 = table2array(test4);
figure
grid on
hold on
yyaxis left
plot(test4(:,3),test4(:,2))
xlabel('Time (Data point taken)') 
ylabel('Current A') 
title('High Speed Test with Default Motor RS 540')
yyaxis right
plot(test4(:,3),test4(:,1))


%% 1A

figure
grid on
hold on
yyaxis left
ylabel('Current A') 
plot(test4(:,3),test4(:,2))
xlabel('Time (Data point taken)')  
yyaxis right
ylabel('IR Distance') 
title('High Speed Test (A) Shows Extreme Currents')
plot(test4(:,3),test4(:,1))
xlim([3400 3600])
% point of interest between 3502-3519

Current1A = mean(test4(3502:3519,2))
Speed1A = Distance/((3519-3502)*SampleSpeed)

%% 1B
figure
grid on
hold on
yyaxis left
ylabel('Current A') 
plot(test4(:,3),test4(:,2))
xlabel('Time (Data point taken)')  
yyaxis right
ylabel('IR Distance') 
title('High Speed Test (B) Inconsistency In Current and IR Readings')

plot(test4(:,3),test4(:,1))
xlim([1850 2050])
% speed point of interest between 1964-1984
% speed point of interest between 1956-1972

Current1B = mean(test4(1956:1972,2))
Speed1B = Distance/((1984-1964)*SampleSpeed)

% %% 1C
% figure
% grid on
% hold on
% yyaxis left
% ylabel('Current A') 
% plot(test4(:,3),test4(:,2))
% xlabel('Time (Data point taken)')  
% yyaxis right
% ylabel('IR Distance') 
% title(['Test 1A'])
% plot(test4(:,3),test4(:,1))
% xlim([950 1150])
% % speed point of interest between 1061-1084
% % speed point of interest between 1956-1972
% 
% Current1C = mean(test4(1061:1084,2))
% Speed1C = Distance/((1084-1061)*SampleSpeed)
%% Official Tests
Current1 = [];
Speed1 = [];
Current2 = [];
Speed2 = [];
Current3 = [];
Speed3 = [];
Current4 = [];
Speed4 = [];
%% test 2 on 1st motor 

SampleSpeed = 0.02025;
Distance = 1.45;
 load('defaultSlow.mat')
Data = table2array(defaultSlow);
figure
grid on
hold on
yyaxis left
plot(Data(:,3),Data(:,2))
xlabel('Time (Data point taken)') 
ylabel('Current A') 
title('Low Speed Test on RS 540 - Default Motor')

yyaxis right
plot(Data(:,3),Data(:,1))

%% 1


CurrentDefSlwA = mean(Data(447:577,2))
SpeedDefSlwA = Distance/((577-447)*SampleSpeed)
EfficiencyDefSlwA = SpeedDefSlwA/CurrentDefSlwA
Current1 =[Current1, CurrentDefSlwA];
Speed1 =[Speed1, SpeedDefSlwA];


%% 2

CurrentDefSlwB = mean(Data(1669:1782,2))
SpeedDefSlwB = Distance/((1785-1669)*SampleSpeed)
EfficiencyDefSlwB = SpeedDefSlwB/CurrentDefSlwB
Current1 =[Current1, CurrentDefSlwB];
Speed1 =[Speed1, SpeedDefSlwB];


%% 3


CurrentDefSlwC = mean(Data(2696:2757,2))
SpeedDefSlwC = Distance/((2757-2696)*SampleSpeed)
EfficiencyDefSlwC = SpeedDefSlwC/CurrentDefSlwC
Current1 =[Current1, CurrentDefSlwC];
Speed1 =[Speed1, SpeedDefSlwC];


%% importing fast simulation
load('defaultFast.mat')
Data = table2array(defaultFast);
figure
grid on
hold on
yyaxis left
plot(Data(:,3),Data(:,2))
xlabel('Time (Data point taken)') 
ylabel('Current A') 
title('High Speed Test (More Controlled) on RS 540 - Default Motor')
yyaxis right
plot(Data(:,3),Data(:,1))

%% 4



CurrentDefSlwD = mean(Data(1192:1224,2))
SpeedDefSlwD = Distance/((1224-1192)*SampleSpeed)
EfficiencyDefSlwD = SpeedDefSlwA/CurrentDefSlwA
Current1 =[Current1, CurrentDefSlwD];
Speed1 =[Speed1, SpeedDefSlwD];


%% 5



CurrentDefSlwE = mean(Data(2359:2382,2))
SpeedDefSlwE = Distance/((2382-2359)*SampleSpeed)
EfficiencyDefSlwE = SpeedDefSlwA/CurrentDefSlwA
Current1 =[Current1, CurrentDefSlwE];
Speed1 =[Speed1, SpeedDefSlwE];


%% 6


CurrentDefSlwF = mean(Data(3804:3840,2))
SpeedDefSlwF = Distance/((3840-3804)*SampleSpeed)
EfficiencyDefSlwF = SpeedDefSlwA/CurrentDefSlwA
Current1 =[Current1, CurrentDefSlwF];
Speed1 =[Speed1, SpeedDefSlwF];



%% Doing the second motor, 55T

load('T1.mat')
Data = table2array(T1);
figure
grid on
hold on
yyaxis left
plot(Data(:,3),Data(:,2))
xlabel('Time (Data point taken)') 
ylabel('Current A') 
title('Varying Speed Test on 55T Motor')

yyaxis right
plot(Data(:,3),Data(:,1))



%% 1B

Current55TA = mean(Data(220:268,2))
Speed55TA = Distance/((268-220)*SampleSpeed)
Efficiency55TA = Speed55TA/Current55TA
Current2 =[Current2, Current55TA];
Speed2 =[Speed2, Speed55TA];
%% 2B

Current55TB = mean(Data(945:998,2))
Speed55TB = Distance/((998-945)*SampleSpeed)
Efficiency55TB = Speed55TB/Current55TB
Current2 =[Current2, Current55TB];
Speed2 =[Speed2, Speed55TB];

%% 3B

Current55TC = mean(Data(1631:1677,2))
Speed55TC = Distance/((1677-1631)*SampleSpeed)
Efficiency55TC = Speed55TC/Current55TC
Current2 =[Current2, Current55TC];
Speed2 =[Speed2, Speed55TC];

%% 4B

Current55TD = mean(Data(2512:2642,2))
Speed55TD = Distance/((2642-2512)*SampleSpeed)
Efficiency55TD = Speed55TD/Current55TD
Current2 =[Current2, Current55TD];
Speed2 =[Speed2, Speed55TD];
%% 5B

Current55TE = mean(Data(3065:3104,2))
Speed55TE = Distance/((3104-3065)*SampleSpeed)
Efficiency55TE = Speed55TE/Current55TE
Current2 =[Current2, Current55TE];
Speed2 =[Speed2, Speed55TE];

%% 6B

Current55TF = mean(Data(3706:3792,2))
Speed55TF = Distance/((3792-3706)*SampleSpeed)
Efficiency55TF = Speed55TF/Current55TF
Current2 =[Current2, Current55TF];
Speed2 =[Speed2, Speed55TF];



%% Doing the third motor, 21T

load('T21.mat')
Data = table2array(T21);
figure
grid on
hold on
yyaxis left
plot(Data(:,3),Data(:,2))
xlabel('Time (Data point taken)') 
title('Varying Speed Test on 21T Motor')
ylabel('Current A') 
yyaxis right
plot(Data(:,3),Data(:,1))



%% 1C

Current21TA = mean(Data(198:239,2))
Speed21TA = Distance/((239-198)*SampleSpeed)
Efficiency21TA = Speed21TA/Current21TA
Current3 =[Current3, Current21TA];
Speed3 =[Speed3, Speed21TA];
%% 2C

Current21TB = mean(Data(1570:1609,2))
Speed21TB = Distance/((1609-1570)*SampleSpeed)
Efficiency21TB = Speed21TB/Current21TB
Current3 =[Current3, Current21TB];
Speed3 =[Speed3, Speed21TB];

%% 3C

Current21TC = mean(Data(2861:2940,2))
Speed21TC = Distance/((2940-2861)*SampleSpeed)
Efficiency21TC = Speed21TC/Current21TC
Current3 =[Current3, Current21TC];
Speed3 =[Speed3, Speed21TC];

%% 4C

Current21TD = mean(Data(4042:4146,2))
Speed21TD = Distance/((4146-4042)*SampleSpeed)
Efficiency21TD = Speed21TD/Current21TD
Current3 =[Current3, Current21TD];
Speed3 =[Speed3, Speed21TD];
%% 5C

Current21TE = mean(Data(5013:5224,2))
Speed21TE = Distance/((5224-5013)*SampleSpeed)
Efficiency21TE = Speed21TE/Current21TE;
Current3 =[Current3, Current21TE];
Speed3 =[Speed3, Speed21TE];

%% 6C

Current21TF = mean(Data(5775:5814,2))
Speed21TF = Distance/((5814-5775)*SampleSpeed)
Efficiency21TF = Speed21TF/Current21TF
Current3 =[Current3, Current21TF];
Speed3 =[Speed3, Speed21TF];

%% Doing the fourth motor, 80T

load('T80.mat')
Data = table2array(T80);
figure
grid on
hold on
yyaxis left
plot(Data(:,3),Data(:,2))
xlabel('Time (Data point taken)') 
title('Varying Speed Test on 80T Motor (A)')

ylabel('Current A') 
yyaxis right
plot(Data(:,3),Data(:,1))



%% 1D

Current80TA = mean(Data(327:450,2))
Speed80TA = Distance/((450-327)*SampleSpeed)
Efficiency80TA = Speed80TA/Current80TA
Current4 =[Current4, Current80TA];
Speed4 =[Speed4, Speed80TA];
%% 2D

Current80TB = mean(Data(1396:1569,2))
Speed80TB = Distance/((1569-1396)*SampleSpeed)
Efficiency80TB = Speed80TB/Current80TB
Current4 =[Current4, Current80TB];
Speed4 =[Speed4, Speed80TB];

%% 3D

Current80TC = mean(Data(2453:2511,2))
Speed80TC = Distance/((2511-2453)*SampleSpeed)
Efficiency80TC = Speed80TC/Current80TC
Current4 =[Current4, Current80TC];
Speed4 =[Speed4, Speed80TC];

%% 4D

Current80TD = mean(Data(3319:3369,2))
Speed80TD = Distance/((3369-3319)*SampleSpeed)
Efficiency80TD = Speed80TD/Current80TD
Current4 =[Current4, Current80TD];
Speed4 =[Speed4, Speed80TD];
%% 5D

Current80TE = mean(Data(4212:4302,2))
Speed80TE = Distance/((4302-4212)*SampleSpeed)
Efficiency80TE = Speed80TE/Current80TE
Current4 =[Current4, Current80TE];
Speed4 =[Speed4, Speed80TE];



%% Doing the fourth motor, 80TB

load('T80b.mat')
Data = table2array(T80b);
figure
grid on
hold on
yyaxis left
plot(Data(:,3),Data(:,2))
xlabel('Time (Data point taken)') 
title('Varying Speed Test on 80T Motor  (B)')

ylabel('Current A') 
yyaxis right
plot(Data(:,3),Data(:,1))

%% 6D

Current80TF = mean(Data(219:269,2))
Speed80TF = Distance/((269-219)*SampleSpeed)
Efficiency80TF = Speed80TF/Current80TF
Current4 =[Current4, Current80TF];
Speed4 =[Speed4, Speed80TF];

%% 7D

Current80TG = mean(Data(1337:1382,2))
Speed80TG = Distance/((1382-1337)*SampleSpeed)
Efficiency80TF = Speed80TG/Current80TG
Current4 =[Current4, Current80TG];
Speed4 =[Speed4, Speed80TG];
%% 8D

Current80TH = mean(Data(2315:2353,2))
Speed80TH = Distance/((2353-2315)*SampleSpeed)
Efficiency80TH = Speed80TH/Current80TH
Current4 =[Current4, Current80TH];
Speed4 =[Speed4, Speed80TH];
%% 9D

Current80TH = mean(Data(3216:3266,2))
Speed80TH = Distance/((3266-3216)*SampleSpeed)
Efficiency80TH = Speed80TF/Current80TH
Current4 =[Current4, Current80TH];
Speed4 =[Speed4, Speed80TH];




%%
figure 
hold on
scatter(Current1',Speed1','b')
scatter(Current2',Speed2','c')
scatter(Current3',Speed3','r')
scatter(Current4',Speed4','g')

b1 = Current1'\Speed1';
plot(Current1',Current1'*b1,'b')
b2 = Current2'\Speed2';
plot(Current2',Current2'*b2,'c')
b3 = Current3'\Speed3';
plot(Current3',Current3'*b3,'r')
b4 = Current4'\Speed4';
plot(Current4',Current4'*b4,'g')





xlabel('Average Current A')  
ylabel('Speed m/s')  
xlim([0 3.5])
ylim([0 3])
legend('540/90W','540/55T','540/21T','540/80T','540/90W Reg','540/55T Reg','540/21T Reg','540/80T Reg')
title('Comparing Efficiency of Varying Motors for Varying Induced Current')
Current1 = [Current1,0,0,0];
Current2 = [Current2,0,0,0];
Current3 = [Current3,0,0,0];
Speed1 = [Speed1,0,0,0];
Speed2 = [Speed2,0,0,0];
Speed3 = [Speed3,0,0,0];



%TestCurrent = [Current1',Current2',Current3',Current4'];
%TestSpeed = [Speed1',Speed2',Speed3',Speed4'];

%% testing speeds


load('TEST7.mat')
Data = table2array(TEST7);
figure
grid on
hold on

MovAvg = movavg(Data(50:730,2),'simple',20);

plot(Data(50:730,3),Data(50:730,2))
plot(Data(50:730,3),MovAvg,'LineWidth',3.0)
xlabel('Time (Data point taken)') 
ylabel('Current A')
title('Race Simulation Displaying Current  on Acceleration/Straights/Corners/Stopping')
xlim([50 730])
MeanCurrent = mean(Data(50:730,2))



