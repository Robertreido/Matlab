clear all; 
close all;
clc;
clf

tic
%% Project 1: Distributed MIMO
%% Up Link MRC B_k = 1
for fold = 1
M = 6; %BS Antennas
K = 4; %Number of uses
B_k = 1;%Power received
nt = 1;
N =30;%N of SNR values
Nb = 1000; %Number of simulations
SNR =  linspace(-20,20,N);%SNR in dB, SNR p_ul/o^2
SNR_linear = db2pow(SNR);
SigmaSquared = 1;%Noise Power
SINRk_avg = zeros([N , K]);
R = zeros([numel(SNR) 4 Nb]);


for iter = 1:numel(SNR)
    p_ul = SigmaSquared*SNR_linear(iter);
    for i = 1:Nb
        h = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        for k = 1:K
            Sum = 0;
            for i1 = 1:K
                if i1 ~= k
                    temp =  B_k * abs(det(h(:,k)'*h(:,i1)))^2;
                    Sum = Sum + temp;
                end
            end         
            SINRnumer = (p_ul*B_k*norm(h(:,k))^4);
            SINRdenom = (SigmaSquared*norm(h(:,k))^2 +p_ul*Sum);
            SINRk = SINRnumer/SINRdenom;
            R(iter,k,i) = log2(1+SINRk);
        end
    end
end
Rmean = mean(R,3);
Rsum = sum(Rmean,2);

figure(1)
hold on
grid on
plot(SNR,Rsum,'-.r','LineWidth',2)

"MRC UL done (1)"
% ****************************************************
% ZF Combining

R = zeros([numel(SNR) 4 Nb]);
for iter = 1:numel(SNR)
    p_ul= SigmaSquared*SNR_linear(iter);
    for i = 1:Nb
        h = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        for k = 1:K
            temp = inv(h'*h);
            SINRk = (p_ul*B_k)/(SigmaSquared*temp(k,k));
            R(iter,k,i) = log2(1+SINRk);
        end
        
    end
end
Rmean = mean(R,3);
Rsum = sum(Rmean,2);


plot(SNR,Rsum,'--','color',[0 0 0.5],'LineWidth',2)

"ZF UL done (1)"
% ****************************************************
% MMSE

SNR =  linspace(-20,20,N);%SNR in dB, SNR p_ul/o^2
SNR_linear = db2pow(SNR);
SigmaSquared = 1;
SINRk_avg = zeros([N , K]);
R = zeros([numel(SNR) 4 Nb]);
Ident = eye(M);

for iter = 1:numel(SNR)
    p_ul= SigmaSquared*SNR_linear(iter);
    for i = 1:Nb
        h = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        
        for k = 1:K
            Sum = 0;
            
            for i1 = 1:K
                if i1 ~= k
                    temp =  h(:,i1)* h(:,i1)';
                    Sum = Sum + temp;
                end
            end
            
            SINRk = p_ul*B_k*h(:,k)'*inv(p_ul*B_k*Sum + Ident)*h(:,k);
            R(iter,k,i) = log2(1+SINRk);
        end
        
    end
end
Rmean = mean(R,3);
Rsum = abs(sum(Rmean,2));

"MMSE UL done(1)"

plot(SNR,Rsum,'color',[0 0.5 0],'LineWidth',2)
grid on
xlabel('SNR dB');
ylabel('Sum Rate (bits/Hz)');
legend('MRC','ZF','MMSE')
title('Performance of linear receivers in the uplink, B_k = 1')
ylim([0 35])

"Uplink done(1)"

end
%% Downlink
for fold = 1
% MRC
M = 6; %BS Antennas
K = 4; %Number of uses
B_k = 1;
N =30;
Nb = 1000; %averaging N
SNR =  linspace(-10,30,N);%SNR in dB, SNR p_ul/o^2
SNR_linear = db2pow(SNR);
SigmaSquared = 1;


%calculating eta
Trace = [];
for i = 1:10000
    H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H;
        Trace(i,1) = trace(W*W');
end
eta = 1/mean(Trace);



R = zeros([numel(SNR) 4 Nb]);
for iter = 1:numel(SNR)
    p_dl = SigmaSquared*SNR_linear(iter);
    for i = 1:Nb      
        H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H;

         for k = 1:K
            Sum = 0;
            for i1 = 1:K
                if i1 ~= k
                    temp = abs(det(H(:,k)'*W(:,i1)))^2;
                    Sum = Sum + temp;                    
                end
            end
            
            SINRnumer = (p_dl*eta*B_k*(det(H(:,k)'*W(:,k)))^2);
            SINRdenom = (SigmaSquared + p_dl*eta*B_k*Sum);
            SINRk = SINRnumer/SINRdenom;
            R(iter,k,i) = log2(1+SINRk);
        end
    end
end
Rmean = mean(R,3);
Rsum = sum(Rmean,2);

figure(2)
hold on
grid on
plot(SNR,Rsum,'-.r','LineWidth',2)
ylim([0 40])

"MRC DL done (1)"
% ****************************************************
% ZF


%calculating eta
Trace = zeros(1000,1);
for i = 1:10000
    H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
    W = H*inv(H'*H);
    Trace(i,1) = trace(W*W');
end
eta = 1./mean(Trace);
eta_ZF = eta;
eta = 0.35;


R = zeros([numel(SNR) 4 Nb]);
for iter = 1:numel(SNR)
    p_dl = SigmaSquared*SNR_linear(iter);
    for i = 1:Nb      
        H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H*inv(H'*H);

         for k = 1:K
            Sum = 0;
            for i1 = 1:K
                if i1 ~= k
                    temp = (det(H(:,k)'*W(:,i1)))^2;
                    Sum = Sum + temp;                    
                end
            end
            
            SINRnumer = (p_dl*eta*B_k*(det(H(:,k)'*W(:,k)))^2);
            SINRdenom = (SigmaSquared + p_dl*eta*B_k*Sum);
            SINRk = SINRnumer/SINRdenom;
            R(iter,k,i) = log2(1+SINRk);
        end
    end
end
Rmean = mean(R,3);
Rsum = abs(sum(Rmean,2));

hold on
grid on
plot(SNR,Rsum,'--b','LineWidth',2)
ylim([0 40])

"ZF DL done (1)"
% ****************************************************
% MMSE

%calculating eta
Trace = zeros(numel(SNR),N*10);
Ident = eye(K);
    for i = 1:N*10  
        for iter = 1:numel(SNR)
        H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H*inv(H'* H + Ident * K/SNR_linear(iter));
        %W = SNR_linear(iter);
        Trace(iter,i) = trace(W*W');
    end
end 
temp = mean(Trace');
eta = (1./temp);
eta_MMSE = eta;
%figure()
%scatter(SNR,eta)


%


for iter = 1:numel(SNR)
    p_dl = SigmaSquared*SNR_linear(iter);
    for i = 1:Nb      
        H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H* inv(H'*H + Ident*K./SNR_linear(iter));

         for k = 1:K
            Sum = 0;
            for i1 = 1:K
                if i1 ~= k
                    temp = (det(H(:,k)'*W(:,i1)))^2;
                    Sum = Sum + temp;                    
                end
            end
            SINRnumer = (p_dl*eta(iter)*B_k*(det(H(:,k)'*W(:,k)))^2);
            SINRdenom = (SigmaSquared + p_dl*eta(iter)*B_k*Sum);            
            SINRk = SINRnumer/SINRdenom;
            R(iter,k,i) = log2(1+SINRk);
        end
    end
end
Rmean = mean(R,3);
Rsum = abs(sum(Rmean,2));

hold on
grid on
plot(SNR,Rsum,'color',[0 0.5 0],'LineWidth',2)
ylim([0 40])
"MMSE DL done (1)"

xlabel('SNR dB');
ylabel('Sum Rate (bits/Hz)');
legend('MRC','ZF','MMSE')
title('Performance of linear precoders in the downlink, B_k = 1')



end
%% CDF determing d
for fold = 1
 N = 30;
 d = 1.68; 
 gamma = 3;
 Pr = [];
 Pt = 20; % Transmit Power
Pt_pow = db2pow(Pt)'; 
%L =mean(10.^(normrnd(0,16,10000)/10),2);
L =  (10.^(normrnd(0,8,10000,1)/10));
Pr= Pt_pow*L*d^(-gamma);

Pr = (pow2db(Pr));
figure (3)
cdfplot(Pr)
txt = {'Boundary = 1.68 Units'};
text(-10,0.75,txt);
xlabel('Power Received dB');
ylabel('Percentage');
title('CDF to determine 95% prob of Power_r > 0dB to find Boundary')
"CDF determing D done"
end
%% Up Link Varying Beta_k
clear
for fold = 1
%MRC
M = 6; %BS Antennas
K = 4; %Number of uses
B_k = 1;
nt = 1;
N =30;
Nb = 10000; %averaging N
SNR =  linspace(-20,20,N);%SNR in dB, SNR p_ul/o^2
SigmaSquared = 1;
R = zeros([numel(SNR) 4 Nb]);
d = 1.68;; 
gamma = 3; %pathloss exponent
SNRpow = db2pow(SNR);

for iter = 1:numel(SNR)
    p_ul = SigmaSquared*SNRpow(iter);
    for i = 1:Nb
        h = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
        Phases = -d +(2*d)*rand(K,1);
        d_user = (Mags .* exp(j*Phases));
        distance= (abs(d_user)*[1 1 1 1 1 1])';
        L = (10.^(normrnd(0,8)/10));
        Beta = SNRpow(iter)*L.*distance.^(-gamma);
        
        for k = 1:K
            Sum = 0;
            for i1 = 1:K
                if i1 ~= k
                    temp =  Beta(:,k) * abs(det(h(:,k)'*h(:,i1)))^2;
                    Sum = Sum + temp;
                end
            end
            
            SINRnumer = (p_ul*Beta(:,k)*norm(h(:,k))^4);
            SINRdenom = (SigmaSquared*norm(h(:,k))^2 +p_ul*Sum);
            SINRk = SINRnumer./SINRdenom;
            R(iter,k,i) = log2(1+mean(SINRk));
        end
    end
end
Rmean = mean(R,3);
Rsum = sum(Rmean,2);

figure(4)
hold on
grid on
plot(SNR,Rsum,'-.r','LineWidth',2)
ylim([0 40])
"MRC UL done (2)"
% ****************************************************
%% ZF Combining

R = zeros([numel(SNR) 4 Nb]);
for iter = 1:numel(SNR)
    p_ul= SigmaSquared*SNRpow(iter);
    for i = 1:Nb
        h = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
         Phases = -d +(2*d)*rand(K,1);
         d_user = (Mags .* exp(j*Phases));
         distance= (abs(d_user)*[1 1 1 1 1 1])';
         L = (10.^(normrnd(0,8)/10));         
         Beta = SNRpow(iter)*L.*distance.^(-gamma);
        for k = 1:K
            temp = inv(h'*h);
            SINRk = (p_ul*Beta(:,k))./(SigmaSquared*temp(k,k));
            R(iter,k,i) = log2(1+mean(SINRk));
        end
        
    end
end
Rmean = mean(R,3);
Rsum = sum(Rmean,2);


plot(SNR,Rsum,'--','color',[0 0 0.5],'LineWidth',2)

"ZF UL done (2)"
% ****************************************************
% MMSE

R = zeros([numel(SNR) 4 Nb]);
Ident = eye(M);

for iter = 1:numel(SNR)
    p_ul= SigmaSquared*SNRpow(iter);
    for i = 1:Nb
        h = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
        Phases = -d +(2*d)*rand(K,1);
        d_user = (Mags .* exp(j*Phases));
        distance= (abs(d_user)*[1 1 1 1 1 1])';
        L = (10.^(normrnd(0,8)/10));
        Beta = SNRpow(iter)*L.*distance.^(-gamma);
        for k = 1:K
            Sum = 0;     
            for i1 = 1:K
                if i1 ~= k
                    temp =  h(:,i1)* h(:,i1)';
                    Sum = Sum + temp;
                end
            end
            SINRk = p_ul*Beta(:,k)*h(:,k)'*inv(p_ul*Beta(:,k).*Sum + Ident)*h(:,k);
            R(iter,k,i) = log2(1+mean(SINRk));
        end
        
    end
end
Rmean = mean(R,3);
Rsum = abs(sum(Rmean,2));

"MMSE UL done (2)"

plot(SNR,Rsum,'color',[0 0.5 0],'LineWidth',2)
grid on
xlabel('SNR dB');
ylabel('Sum Rate (bits/Hz)');
legend('MRC','ZF','MMSE')
title('Performance of linear receivers in the uplink, Varying B_k')
%ylim([0 35])

"Uplink done (2)"

end
%% Downlink Varying B_k
clear
for fold = 1
% MRC
M = 6; %BS Antennas
K = 4; %Number of uses
B_k = 1;
N =30;
Nb = 5000; %averaging N
SNR =  linspace(-10,20,N);%SNR in dB, SNR p_ul/o^2
SNR_linear = db2pow(SNR);
SigmaSquared = 1;

 d = 1.68; 
 gamma = 3; %pathloss exponent

%calculating eta
Trace = [];
for i = 1:10000
    H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H;
        Trace(i,1) = trace(W*W');
end
eta = 1/mean(Trace);
SNRpow= db2pow(SNR);

R = zeros([numel(SNR) 4 Nb]);
for iter = 1:numel(SNR)
    p_dl = SigmaSquared*SNR_linear(iter);
    for i = 1:Nb      
        H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H;
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
         Phases = -d +(2*d)*rand(K,1);
         d_user = (Mags .* exp(j*Phases));
         distance= (abs(d_user)*[1 1 1 1 1 1])';
         L = (10.^(normrnd(0,8)/10));         
         Beta = SNRpow(iter)*L.*distance.^(-gamma);
         
         for k = 1:K
            Sum = 0;
            for i1 = 1:K
                if i1 ~= k
                    temp = abs(det(H(:,k)'*W(:,i1)))^2;
                    Sum = Sum + temp;                    
                end
            end
            
            SINRnumer = (p_dl*eta*Beta(:,k)*(det(H(:,k)'*W(:,k)))^2);
            SINRdenom = (SigmaSquared + p_dl*eta*Beta(:,k)*Sum);
            SINRk = mean(SINRnumer./SINRdenom);
            R(iter,k,i) = log2(1+SINRk);
        end
    end
end
Rmean = mean(R,3);
Rsum = sum(Rmean,2);

figure(5)
hold on
grid on
plot(SNR,Rsum,'-.r','LineWidth',2)
%ylim([0 40])

"MRC DL done (2)"
% ****************************************************
%% ZF


%calculating eta
Trace = zeros(1000,1);
for i = 1:10000
    H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
    W = H*inv(H'*H);
    Trace(i,1) = trace(W*W');
end
eta = 1./mean(Trace);
eta_ZF = eta;
eta = 0.35;

R = zeros([numel(SNR) 4 Nb]);
for iter = 1:numel(SNR)
    p_dl = SigmaSquared*SNR_linear(iter);
    for i = 1:Nb      
        H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H*inv(H'*H);
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
         Phases = -d +(2*d)*rand(K,1);
         d_user = (Mags .* exp(j*Phases));
         distance= (abs(d_user)*[1 1 1 1 1 1])';
         L = (10.^(normrnd(0,8)/10));         
         Beta = SNRpow(iter)*L.*distance.^(-gamma);

         for k = 1:K
            Sum = 0;
            for i1 = 1:K
                if i1 ~= k
                    temp = (det(H(:,k)'*W(:,i1)))^2;
                    Sum = Sum + temp;                    
                end
            end
            
            SINRnumer = (p_dl*eta*Beta(:,k)*(det(H(:,k)'*W(:,k)))^2);
            SINRdenom = (SigmaSquared + p_dl*eta*Beta(:,k).*Sum);
            SINRk = mean(SINRnumer./SINRdenom);
            R(iter,k,i) = log2(1+SINRk);
        end
    end
end
Rmean = mean(R,3);
Rsum = abs(sum(Rmean,2));

hold on
grid on
plot(SNR,Rsum,'--b','LineWidth',2)
%ylim([0 40])

"ZF DL done (2)"
% ****************************************************
% MMSE

%calculating eta
Trace = zeros(numel(SNR),N*10);
Ident = eye(K);
    for i = 1:N*10  
        for iter = 1:numel(SNR)
        H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H*inv(H'* H + Ident * K/SNR_linear(iter));
        %W = SNR_linear(iter);
        Trace(iter,i) = trace(W*W');
    end
end 
temp = mean(Trace');
eta = (1./temp);
eta_MMSE = eta;
R = zeros([numel(SNR) 4 Nb]);
for iter = 1:numel(SNR)
    p_dl = SigmaSquared*SNR_linear(iter);
    for i = 1:Nb      
        H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H* inv(H'*H + Ident*K./SNR_linear(iter));
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
         Phases = -d +(2*d)*rand(K,1);
         d_user = (Mags .* exp(j*Phases));
         distance= (abs(d_user)*[1 1 1 1 1 1])';
         L = (10.^(normrnd(0,8)/10));         
         Beta = SNRpow(iter)*L.*distance.^(-gamma);

         for k = 1:K
            Sum = 0;
            for i1 = 1:K
                if i1 ~= k
                    temp = (det(H(:,k)'*W(:,i1)))^2;
                    Sum = Sum + temp;                    
                end
            end
            SINRnumer = (p_dl*eta(iter)*Beta(:,k)*(det(H(:,k)'*W(:,k)))^2);
            SINRdenom = (SigmaSquared + p_dl*eta(iter)*Beta(:,k).*Sum);            
            SINRk = mean(SINRnumer./SINRdenom);
            R(iter,k,i) = log2(1+SINRk);
        end
    end
end
Rmean = mean(R,3);
Rsum = abs(sum(Rmean,2));

hold on
grid on
plot(SNR,Rsum,'color',[0 0.5 0],'LineWidth',2)
%ylim([0 40])
"MMSE DL done (2)"

xlabel('SNR dB');
ylabel('Sum Rate (bits/Hz)');
legend('MRC','ZF','MMSE')
title('Performance of linear precoders in the downlink, Varying B_k')


end
%% Geographical Location Figures
clear
for fold = 1
N =4;
R =1.68;
figure
hold on
scatter([0],[0],100,'^','LineWidth',2)
ang=0:0.01:2*pi; 
xp=R*cos(ang);
yp=R*sin(ang);
plot(xp,yp,'--');
for inc = 1:10
Mags = sqrt(-R^2 +(2*R^2)*rand(N,1));
Phases = -R +(2*R)*rand(N,1);
d = (Mags .* exp(j*Phases));
scatter(imag(d), real(d),'+','LineWidth',2)
end
grid on
xlabel('X units');
ylabel('Y units');
legend('Base Station','Cell Boundary, 1.68 units','Random User Groups (of 4)')
title('System with a Singular Basestation')


figure
hold on
scatter([-R/sqrt(2),R/sqrt(2)],[0,0],100,'^','LineWidth',2)
grid on
ang=0:0.01:2*pi; 
xp=R*cos(ang);
yp=R*sin(ang);
plot(xp,yp,'--');
for inc = 1:10
Mags = sqrt(-R^2 +(2*R^2)*rand(N,1));
Phases = -R +(2*R)*rand(N,1);
d = (Mags .* exp(j*Phases));
scatter(imag(d), real(d),'+','LineWidth',2)
end
xlabel('X units');
ylabel('Y units');
legend('Base Stations','Cell Boundary, 1.68 units','Random User Groups (of 4)')
title('System with Two Basestations')


figure
hold on
M=3;
Temp = R/sqrt(2)*cos(linspace(0,2*pi,M+1))+i*R/sqrt(2)*sin(linspace(0,2*pi,M+1));
Temp = Temp(1:M);
scatter(real(Temp),imag(Temp),100,'^','LineWidth',2)
grid on
ang=0:0.01:2*pi; 
xp=R*cos(ang);
yp=R*sin(ang);
plot(xp,yp,'--');
for inc = 1:10
Mags = sqrt(-R^2 +(2*R^2)*rand(N,1));
Phases = -R +(2*R)*rand(N,1);
d = (Mags .* exp(j*Phases));
scatter(imag(d), real(d),'+','LineWidth',2)
end
xlabel('X units');
ylabel('Y units');
legend('Base Stations','Cell Boundary, 1.68 units','Random User Groups (of 4)')
title('System with  Three Basestations')

"Geographical Diagrams Done"
end

%% Multiple Base Stations
%% Uplink with Multiple BS
clear
for fold = 1
%MRC
M = 6; %BS Antennas
K = 4; %Number of uses
B_k = 1;
nt = 1;
N =30;
Nb = 10000; %averaging N
SNR =  linspace(-20,20,N);%SNR in dB, SNR p_ul/o^2
SNR_linear = db2pow(SNR);
SigmaSquared = 1;
%p_ul = 1;
SINRk_avg = zeros([N , K]);
R = zeros([numel(SNR) 4 Nb]);
count = 0;

 d = 1.68;%33; 
 gamma = 3; %pathloss exponent
SNRpow = db2pow(SNR);

Temp = d/sqrt(2)*cos(linspace(0,2*pi,M+1))+i*d/sqrt(2)*sin(linspace(0,2*pi,M+1));
Bs_loc = Temp(1:M);
for iter = 1:numel(SNR)
    p_ul = SigmaSquared*SNR_linear(iter);
    for inc = 1:Nb
        % SNR_linear = db2pow(SNR(1,iter));
        h = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
         Phases = -d +(2*d)*rand(K,1);
         d_user = (Mags .* exp(j*Phases));
         distance = abs(d_user-Bs_loc)';
         L = (10.^(normrnd(0,8)/10));         
         Beta = SNRpow(iter)*L.*distance.^(-gamma);
         
        for k = 1:K
            Sum = 0;
            for i1 = 1:K
                if i1 ~= k
                    temp =  Beta(:,k) * abs(det(h(:,k)'*h(:,i1)))^2;
                    Sum = Sum + temp;                    
                end
            end
            
            SINRnumer = (p_ul*Beta(:,k)*norm(h(:,k))^4);
            SINRdenom = (SigmaSquared*norm(h(:,k))^2 +p_ul*Sum);
            SINRk = SINRnumer./SINRdenom;
            R(iter,k,inc) = log2(1+mean(SINRk));
        end
    end
end
Rmean = mean(R,3);
Rsum = sum(Rmean,2);

figure(4)
hold on
grid on
plot(SNR,Rsum,'-.r','LineWidth',2)

"MRC UL done (3)"
% ****************************************************
%% ZF Combining

R = zeros([numel(SNR) 4 Nb]);
Temp = d/sqrt(2)*cos(linspace(0,2*pi,M+1))+i*d/sqrt(2)*sin(linspace(0,2*pi,M+1));
Bs_loc = Temp(1:M);
for iter = 1:numel(SNR)
    p_ul= SigmaSquared*SNR_linear(iter);
    for inc = 1:Nb
        h = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
         Phases = -d +(2*d)*rand(K,1);
         d_user = (Mags .* exp(j*Phases));
         distance = abs(d_user-Bs_loc)';
         L = (10.^(normrnd(0,8)/10));         
         Beta = SNRpow(iter)*L.*distance.^(-gamma);
        for k = 1:K
            temp = inv(h'*h);
            SINRk = (p_ul*Beta(:,k))./(SigmaSquared*temp(k,k));
            R(iter,k,inc) = log2(1+mean(SINRk));
        end
        
    end
end
Rmean = mean(R,3);
Rsum = sum(Rmean,2);


plot(SNR,Rsum,'--','color',[0 0 0.5],'LineWidth',2)

"ZF UL done (3)"
% ****************************************************
% MMSE

SINRk_avg = zeros([N , K]);
R = zeros([numel(SNR) 4 Nb]);
Ident = eye(M);
Temp = d/sqrt(2)*cos(linspace(0,2*pi,M+1))+i*d/sqrt(2)*sin(linspace(0,2*pi,M+1));
Bs_loc = Temp(1:M);
for iter = 1:numel(SNR)
    p_ul= SigmaSquared*SNR_linear(iter);
    for inc = 1:Nb
        h = 1/sqrt(2)*(randn(M,K)+i*randn(M,K));
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
         Phases = -d +(2*d)*rand(K,1);
         d_user = (Mags .* exp(j*Phases));
         distance = abs(d_user-Bs_loc)';
         L = (10.^(normrnd(0,8)/10));         
         Beta = SNRpow(iter)*L.*distance.^(-gamma);
        
        for k = 1:K
            Sum = 0;
            
            for i1 = 1:K
                if i1 ~= k
                    temp =  h(:,i1)* h(:,i1)';
                    Sum = Sum + temp;
                end
            end
            
            SINRk = p_ul*Beta(:,k)*h(:,k)'*inv(p_ul*Beta(:,k).*Sum + Ident)*h(:,k);
            R(iter,k,inc) = log2(1+mean(SINRk));
        end
        
    end
end
Rmean = mean(R,3);
Rsum = abs(sum(Rmean,2));

"MMSE UL done (3)"

plot(SNR,Rsum,'color',[0 0.5 0],'LineWidth',2)
grid on
xlabel('SNR dB');
ylabel('Sum Rate (bits/Hz)');
legend('MRC','ZF','MMSE')
title('Uplink Performance with 3 Distributed Basestations')
%ylim([0 35])

"Uplink done (3)"

end

%% Downlink with Multiple BS
clear
for fold = 1
% MRC
M = 6; %BS Antennas
K = 4; %Number of uses
B_k = 1;
N =30;
Nb = 5000; %averaging N
SNR =  linspace(-10,20,N);%SNR in dB, SNR p_ul/o^2
SNR_linear = db2pow(SNR);
SigmaSquared = 1;

  %SNR = linspace(-20,20,N);
 d = 1.68; 
 gamma = 3; %pathloss exponent

%calculating eta
Trace = [];
for inc = 1:10000
    H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H;
        Trace(inc,1) = trace(W*W');
end
eta = 1/mean(Trace);
SNRpow= db2pow(SNR);


Temp = d/sqrt(2)*cos(linspace(0,2*pi,M+1))+i*d/sqrt(2)*sin(linspace(0,2*pi,M+1));
Bs_loc = Temp(1:M);
R = zeros([numel(SNR) 4 Nb]);
for iter = 1:numel(SNR)
    p_dl = SigmaSquared*SNR_linear(iter);
    for inc = 1:Nb      
        H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H;
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
         Phases = -d +(2*d)*rand(K,1);
         d_user = (Mags .* exp(j*Phases));
         distance = abs(d_user-Bs_loc)';
         L = (10.^(normrnd(0,8)/10));         
         Beta = SNRpow(iter)*L.*distance.^(-gamma);
    
         for k = 1:K
            Sum = 0;
            for i1 = 1:K
                if i1 ~= k
                    temp = abs(det(H(:,k)'*W(:,i1)))^2;
                    Sum = Sum + temp;                    
                end
            end
            
            SINRnumer = (p_dl*eta*Beta(:,k)*(det(H(:,k)'*W(:,k)))^2);
            SINRdenom = (SigmaSquared + p_dl*eta*Beta(:,k)*Sum);
            SINRk = mean(SINRnumer./SINRdenom);
            R(iter,k,inc) = log2(1+SINRk);
        end
    end
end
Rmean = mean(R,3);
Rsum = sum(Rmean,2);

figure
hold on
grid on
plot(SNR,Rsum,'-.r','LineWidth',2)
%ylim([0 40])

"MRC DL done (3)"
% ****************************************************
%% ZF


%calculating eta
Trace = zeros(1000,1);
for inc = 1:10000
    H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
    W = H*inv(H'*H);
    Trace(inc,1) = trace(W*W');
end
eta = 1./mean(Trace);
eta_ZF = eta;
eta = 0.35;

Temp = d/sqrt(2)*cos(linspace(0,2*pi,M+1))+i*d/sqrt(2)*sin(linspace(0,2*pi,M+1));
Bs_loc = Temp(1:M);
R = zeros([numel(SNR) 4 Nb]);
for iter = 1:numel(SNR)
    p_dl = SigmaSquared*SNR_linear(iter);
    for inc = 1:Nb      
        H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H*inv(H'*H);
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
         Phases = -d +(2*d)*rand(K,1);
         d_user = (Mags .* exp(j*Phases));
         distance = abs(d_user-Bs_loc)';
         L = (10.^(normrnd(0,8)/10));         
         Beta = SNRpow(iter)*L.*distance.^(-gamma);

         for k = 1:K
            Sum = 0;
            for i1 = 1:K
                if i1 ~= k
                    temp = (det(H(:,k)'*W(:,i1)))^2;
                    Sum = Sum + temp;                    
                end
            end
            
            SINRnumer = (p_dl*eta*Beta(:,k)*(det(H(:,k)'*W(:,k)))^2);
            SINRdenom = (SigmaSquared + p_dl*eta*Beta(:,k).*Sum);
            SINRk = mean(SINRnumer./SINRdenom);
            R(iter,k,inc) = log2(1+SINRk);
        end
    end
end
Rmean = mean(R,3);
Rsum = abs(sum(Rmean,2));

hold on
grid on
plot(SNR,Rsum,'--b','LineWidth',2)
%ylim([0 40])

"ZF DL done (3)"
% ****************************************************
% MMSE

%calculating eta
Trace = zeros(numel(SNR),N*10);
Ident = eye(K);
for inc = 1:N*10
    for iter = 1:numel(SNR)
        H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H*inv(H'* H + Ident * K/SNR_linear(iter));
        %W = SNR_linear(iter);
        Trace(iter,inc) = trace(W*W');
    end
end
temp = mean(Trace');
eta = (1./temp);
eta_MMSE = eta;

Temp = d/sqrt(2)*cos(linspace(0,2*pi,M+1))+i*d/sqrt(2)*sin(linspace(0,2*pi,M+1));
Bs_loc = Temp(1:M);
R = zeros([numel(SNR) 4 Nb]);
for iter = 1:numel(SNR)
    p_dl = SigmaSquared*SNR_linear(iter);
    for inc = 1:Nb      
        H = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        W = H* inv(H'*H + Ident*K./SNR_linear(iter));
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
         Phases = -d +(2*d)*rand(K,1);
         d_user = (Mags .* exp(j*Phases));
         distance = abs(d_user-Bs_loc)';
         L = (10.^(normrnd(0,8)/10));         
         Beta = SNRpow(iter)*L.*distance.^(-gamma);

         for k = 1:K
            Sum = 0;
            for i1 = 1:K
                if i1 ~= k
                    temp = (det(H(:,k)'*W(:,i1)))^2;
                    Sum = Sum + temp;                    
                end
            end
            SINRnumer = (p_dl*eta(iter)*Beta(:,k)*(det(H(:,k)'*W(:,k)))^2);
            SINRdenom = (SigmaSquared + p_dl*eta(iter)*Beta(:,k).*Sum);            
            SINRk = mean(SINRnumer./SINRdenom);
            R(iter,k,inc) = log2(1+SINRk);
        end
    end
end
Rmean = mean(R,3);
Rsum = abs(sum(Rmean,2));

hold on
grid on
plot(SNR,Rsum,'color',[0 0.5 0],'LineWidth',2)
%ylim([0 40])
"MMSE DL done (3)"

xlabel('SNR dB');
ylabel('Sum Rate (bits/Hz)');
legend('MRC','ZF','MMSE')
title('Performance of Downlink with 3 Basestations')


end

%% CDF Uplink plots  for MMSE @ SNR = 10dB, M = 6
for fold = 1

M = 6; %BS Antennas
K = 4; %Number of uses
N =30;
Nb = 10000; %averaging N
SNR =  linspace(-20,20,N);%SNR in dB, SNR p_ul/o^2
SNR_linear = db2pow(SNR);
SigmaSquared = 1;%noise Power
SINRk_avg = zeros([N , K]);
R = [];%zeros([1 N Nb]);
count = 0;
d = 1.68;%33; 
gamma = 3; %pathloss exponent
SNRpow = db2pow(SNR);


% MMSE

SNR = 10; %SNR at 10dB is the anaylsis point
SNRpow = db2pow(SNR);
R = zeros([numel(SNR) 4 Nb]);
Ident = eye(M);

figure(25)
hold on

Bs_loc = [];
Bs_loc(1,:) = [0,0,0,0,0,0];
Bs_loc(2,:) = [-i*d/sqrt(2),-i*d/sqrt(2),-i*d/sqrt(2),i*d/sqrt(2),i*d/sqrt(2),i*d/sqrt(2)];
%Bs_loc(3,:) = [-i*d/2,-i*d/2,-i*d/2,i*d/2,i*d/2,i*d/2]/2;
Temp = d/sqrt(2)*cos(linspace(0,2*pi,M+1))+i*d/sqrt(2)*sin(linspace(0,2*pi,M+1));
Bs_loc(3,:) =  Temp(1:M);

for iterb = 1:3
for iter = 1:1
    p_ul= SigmaSquared*SNR_linear(iter);
    for inc = 1:Nb
        h = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));       
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
         Phases = -d +(2*d)*rand(K,1);
         d_user = (Mags .* exp(j*Phases));
         distance = abs(d_user-Bs_loc(iterb,:))';
         L = (10.^(normrnd(0,8)/10));         
         Beta = SNRpow(iter)*L.*distance.^(-gamma);
        
        for k = 1:K
            Sum = 0;
            
            for i1 = 1:K
                if i1 ~= k
                    temp =  h(:,i1)* h(:,i1)';
                    Sum = Sum + temp;
                end
            end
            
            SINRk = p_ul*Beta(:,k)*h(:,k)'*inv(p_ul*Beta(:,k).*Sum + Ident)*h(:,k);
            R(iter,k,inc) =mean(SINRk);% log2(1+mean(SINRk));
        end
        
    end
end
Rcdf = squeeze(sum(abs(R),2));
cdfplot(pow2db(Rcdf))
xlim([-20 40])

end





for iterb = 1:3

for iter = 1:numel(SNR)
    p_ul= SigmaSquared*SNR_linear(iter);
    for inc = 1:Nb
        h = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
         Phases = -d +(2*d)*rand(K,1);
         d_user = (Mags .* exp(j*Phases));
         distance = abs(d_user-Bs_loc(iterb,:))';
         L = (10.^(normrnd(0,8)/10));         
         Beta = SNRpow(iter)*L.*distance.^(-gamma);
        for k = 1:K
            temp = inv(h'*h);
            SINRk = (p_ul*Beta(:,k))./(SigmaSquared*temp(k,k));
            R(iter,k,inc) =mean(SINRk);
        end
        
    end
end
Rcdf = squeeze(sum(abs(R),2));
cdfplot(pow2db(Rcdf))
xlim([-20 40])


end
title('SINR CDF for Varying Basestations for MMSE and ZF')
xlabel('SINR dB')
ylabel('Percentage')
legend('MMSE: BS Origin', 'MMSE: 2 BS','MMSE: 3 BS','ZF: BS Origin','ZF: 2 BS','ZF: 3 BS')

end


%% CDF Uplink plots  for MMSE & ZF @ SNR = 10dB, M = 16
for fold = 1

M = 128; %BS Antennas
K = 4; %Number of uses
N =30;
Nb = 10000; %averaging N
SNR =  linspace(-20,20,N);%SNR in dB, SNR p_ul/o^2
SNR_linear = db2pow(SNR);
SigmaSquared = 1;%noise Power
SINRk_avg = zeros([N , K]);
R = [];%zeros([1 N Nb]);
count = 0;
d = 1.68;%33; 
gamma = 3; %pathloss exponent
SNRpow = db2pow(SNR);


% MMSE

SNR = 10; %SNR at 10dB is the anaylsis point
SNRpow = db2pow(SNR);
R = zeros([numel(SNR) 4 Nb]);
Ident = eye(M);

figure(29)
hold on

d = 1.68;

Bs_loc = [];
Bs_loc(1,:) = zeros([M 1]);
Temp = d/sqrt(2)*cos(linspace(0,2*pi,M+1))+i*d/sqrt(2)*sin(linspace(0,2*pi,M+1));
Temp = Temp(1:M);

Bs_loc(2,:) = Temp;
%Bs_loc(2,:) = [-i*d/2,-i*d/2,-i*d/2,i*d/2,i*d/2,i*d/2];
%Bs_loc(3,:) = [-i*d/2,-i*d/2,-i*d/2,i*d/2,i*d/2,i*d/2]/2;
%Bs_loc(3,:) = [0+i*d/2,0+i*d/2,(-0.433*d)-(i*0.25*d),(-0.433*d)-(i*0.25*d),0.433*d-i*0.25*d,0.433*d-i*0.25*d];

for iterb = 1:2
for iter = 1:1
    p_ul= SigmaSquared*SNR_linear(iter);
    for inc = 1:Nb
        h = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));       
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
         Phases = -d +(2*d)*rand(K,1);
         d_user = (Mags .* exp(j*Phases));
         distance = abs(d_user-Bs_loc(iterb,:))';
         L = (10.^(normrnd(0,8)/10));         
         Beta = SNRpow(iter)*L.*distance.^(-gamma);
        
        for k = 1:K
            Sum = 0;
            
            for i1 = 1:K
                if i1 ~= k
                    temp =  h(:,i1)* h(:,i1)';
                    Sum = Sum + temp;
                end
            end
            
            SINRk = p_ul*Beta(:,k)*h(:,k)'*inv(p_ul*Beta(:,k).*Sum + Ident)*h(:,k);
            R(iter,k,inc) =mean(SINRk);% log2(1+mean(SINRk));
        end
        
    end
end
Rcdf = squeeze(sum(abs(R),2));
cdfplot(pow2db(Rcdf))
xlim([-10 60])

end





for iterb = 1:2

for iter = 1:numel(SNR)
    p_ul= SigmaSquared*SNR_linear(iter);
    for inc = 1:Nb
        h = 1/sqrt(2)*(randn(M,K)+1i*randn(M,K));
        Mags = sqrt(-d^2 +(2*d^2)*rand(K,1));
         Phases = -d +(2*d)*rand(K,1);
         d_user = (Mags .* exp(j*Phases));
         distance = abs(d_user-Bs_loc(iterb,:))';
         L = (10.^(normrnd(0,8)/10));         
         Beta = SNRpow(iter)*L.*distance.^(-gamma);
        for k = 1:K
            temp = inv(h'*h);
            SINRk = (p_ul*Beta(:,k))./(SigmaSquared*temp(k,k));
            R(iter,k,inc) =mean(SINRk);
        end
        
    end
end
Rcdf = squeeze(sum(abs(R),2));
cdfplot(pow2db(Rcdf))
xlim([-10 60])


end

title('SINR CDF for Massive MIMO with 128 Antennas for MMSE and ZF')
xlabel('SINR dB')
ylabel('Percentage')
legend('MMSE: BS Origin', 'MMSE: 128 Distr. BS','ZF: BS Origin','ZF: 128 Distr. BS')
end

toc