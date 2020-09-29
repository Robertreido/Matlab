% Project 2
% Robert Reid

clear; 
close all;
clc;
clf;
%% Determing the Optimal Transmit Power
M = 150; % number of 
C = 3; % number of Cluster
S = 16; % number of subpaths
delta = 0.5; %antenna spacing in wavelengths
Gamma = 3.2; % pathloss exponent
sigma_sf = (8.2); %dB shadow fading standard deviation
sigma_c = 14.4^2; %cluster angle variance
sigma_s = 1.28^2; %subpath angle variance
r = 100;
r_o = 10;
L = 10;
iteration = 10000;
d_o = 1; % reference distance 


p_ul = 70.5;% 70.5dB for 10% at -0dB
p_ul = db2pow(p_ul);
for l = 1:iteration
    Mags =sqrt(r_o^2+ +(r^2-r_o^2)*rand(1,1));
    Phases = -r +(2*r)*rand(1,1);
    d_l = (Mags .* exp(1i*Phases));
    X_l = (10.^(normrnd(0,sigma_sf)/10));
    Beta_l(l) = 1*X_l*(abs(d_l)/d_o).^-Gamma;
    SNR(l) = p_ul * Beta_l(l);
    
        
end
figure(1)
cdfplot(pow2db(SNR))
xlabel('SNR dB');
ylabel('Percentage');
title('CDF calculating probability of SNR > 0dB is > 10\%, Pul = 70.5 dB')
figure(2)
cdfplot(pow2db(Beta_l))
xlabel('Power Received dB');
ylabel('Percentage');
title('CDF of Power Received at Power Transmit = 70.5dB')
%% Displaying Exclusion Zone of Users


N =10;
r =100;
r_o = 10;
figure(10)
hold on
scatter([0],[0],100,'^','LineWidth',2)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(xp,yp,'r-','LineWidth',3);
ang=0:0.01:2*pi; 
xp=r/10*cos(ang);
yp=r/10*sin(ang);
plot(xp,yp,'r-.','LineWidth',3);
for inc = 1:1
    Mags =sqrt(r_o^2+ +(r^2-r_o^2)*rand(N,1));
    Phases = -r +(2*r)*rand(N,1);
    d = (Mags .* exp(1i*Phases));
    scatter(imag(d), real(d),'b+','LineWidth',2)
end
grid on
xlabel('X m');
ylabel('Y m');
legend('Base Station','Cell Boundary, 100 m','Exclusion Boundary 10m','User')
title('Cell with L users')

%% Plotting SINR for vary SNR - Ray Based
SNR = linspace(-30,60,20);
SNRpow = db2pow(SNR);
N = 3000;
L = 10;
M = 150;

p_ul = 70.5;
p_ul = db2pow(p_ul);
A = p_ul; %received power at reference distance
SINR = zeros([N L numel(SNR)]);
Beta_l = zeros([N L]);
phi_c = zeros([1 C]);
Delta_cs = zeros([1 S]);
phi_cs = zeros([C S]);
h_l = zeros([L C S M]);
for snr = 1:numel(SNR)
    for iter = 1:N
        %determine H matrix for L users
        for l = 1:L
            %Dropping users and determining Beta
            Mags =sqrt(r_o^2+ +(r^2-r_o^2)*rand(1,1));
            Phases = -r +(2*r)*rand(1,1);
            d_l = (Mags .* exp(1i*Phases));
            X_l = (10.^(normrnd(0,sigma_sf)/10));% lognormal shadowing
            Beta_l(iter,l) = A*X_l*(abs(d_l)/d_o).^-Gamma;
            phi_c = normrnd(0,sigma_c,C,1);
            
            Betac = Beta_l(iter,l)./(2.^(1:C));
            Betac(C) = Betac(C-1);
            
            for c = 1:C
                Delta_cs = laprnd(S,1,0,sigma_s);
                phi_cs(c,:) = phi_c(c) + Delta_cs;
                Beta_cs_l = Betac(c)/S;
                for s = 1:S               
                    theta_cs_l = (2*pi)*rand(1);
                    gamma_l = sqrt(Beta_cs_l)*exp(i*theta_cs_l);
                    a_l = exp(i*2*pi*(0:M-1)*delta*sind(phi_cs(c,s)))';
                    h_l(l,c,s,:) = gamma_l*a_l;                  
                end
            end
            h_l_sum = sum(h_l,2);
            h_l_sum = squeeze(sum(h_l_sum,3))';                 
        end     
        %determine SINR for each user
        for l = 1:L            
            H = h_l_sum;
            W = H;           
            Sum = 0;
            for lb = 1:L
                if lb ~= l
                    temp = abs(W(:,l)'*H(:,lb))^2;
                    Sum = Sum + temp;
                end
            end

            SINR(iter,l,snr) = (p_ul*(abs(W(:,l)'*H(:,l)))^2)/(p_ul*Sum + (p_ul*Beta_l(iter,l)./SNRpow(snr)) *norm(W(:,l))^2);
        end
    end  
end

SEmrc = squeeze(log2(1+mean(SINR,1)));
SEmrc = sum(SEmrc,1)


SINR = zeros([N L numel(SNR)]);
Beta_l = zeros([N L]);
phi_c = zeros([1 C]);
Delta_cs = zeros([1 S]);
phi_cs = zeros([C S]);
h_l = zeros([L C S M]);
for snr = 1:numel(SNR)
    for iter = 1:N
        %determine H matrix for L users
        for l = 1:L
            Mags =sqrt(r_o^2+ +(r^2-r_o^2)*rand(1,1));
            Phases = -r +(2*r)*rand(1,1);
            d_l = (Mags .* exp(1i*Phases));
            X_l = (10.^(normrnd(0,sigma_sf)/10));
            Beta_l(iter,l) = A*X_l*(abs(d_l)/d_o).^-Gamma;
            phi_c = normrnd(0,sigma_c,C,1);
            
            Betac = Beta_l(iter,l)./(2.^(1:C));
            Betac(C) = Betac(C-1);
            
            for c = 1:C
                Delta_cs = laprnd(S,1,0,sigma_s);
                phi_cs(c,:) = phi_c(c) + Delta_cs;
                Beta_cs_l = Betac(c)/S;
                for s = 1:S                    
                    theta_cs_l = (2*pi)*rand(1);
                    gamma_l = sqrt(Beta_cs_l)*exp(i*theta_cs_l);                
                    a_l = exp(i*2*pi*(0:M-1)*delta*sind(phi_cs(c,s)))';
                    h_l(l,c,s,:) = gamma_l*a_l;                 
                end
            end
            h_l_sum = sum(h_l,2);
            h_l_sum = squeeze(sum(h_l_sum,3))';                   
        end     
        %determine SINR for each user
        for l = 1:L        
            H = h_l_sum;
            W = H*inv(H'*H);
            
            Sum = 0;
            for lb = 1:L
                if lb ~= l
                    temp = abs(W(:,l)'*H(:,lb))^2;
                    Sum = Sum + temp;
                end
            end
            SINR(iter,l,snr) = (p_ul*(abs(W(:,l)'*H(:,l)))^2)/(p_ul*Sum + (p_ul*Beta_l(iter,l)./SNRpow(snr)) *norm(W(:,l))^2);
        end
    end   
end
SEzf = squeeze(log2(1+mean(SINR,1)));
SEzf = sum(SEzf,1);


figure(3)
hold on
plot(SNR,((SEmrc)),'r-')
plot(SNR,(SEzf),'b-')
grid on


%% Kronecker Model
SNR = linspace(-30,60,20)
SNRpow = db2pow(SNR)

N = 3000;

p = 0.9
R = zeros([M L]);
for j = 1:M
    for k = 1:L
        R(j,k) = p^abs(j-k);
    end
end 
SINR = zeros([N L numel(SNR)]);
Beta_l = zeros([N L]);
phi_c = zeros([1 C]);
Delta_cs = zeros([1 S]);
phi_cs = zeros([C S]);
h_l = zeros([L C S M]);
for snr = 1:numel(SNR)
    for iter = 1:N
        %determine H matrix for L users
        for l = 1:L
            
            Mags =sqrt(r_o^2+ +(r^2-r_o^2)*rand(1,1));
            Phases = -r +(2*r)*rand(1,1);
            d_l = (Mags .* exp(1i*Phases));
            X_l = (10.^(normrnd(0,sigma_sf)/10));
            Beta_l(iter,l) = A*X_l*(abs(d_l)/d_o).^-Gamma;
            H = 1/sqrt(2) *(randn(M,L)+1i*randn(M,L));
            H = sqrt(R).*H;
            
        end  
        
        
        %determine SINR for each user
        for l = 1:L
            
            W = H;
            Sum = 0;
            for lb = 1:L
                if lb ~= l
                    temp = abs(W(:,l)'*H(:,lb))^2;
                    Sum = Sum + temp;
                end
            end
            
            SINR(iter,l,snr) = (p_ul*(abs(W(:,l)'*H(:,l)))^2)/(p_ul*Sum + (p_ul*Beta_l(iter,l)./SNRpow(snr)) *norm(W(:,l))^2);
            
        end
    end
    
end

SEmrc = squeeze(log2(1+mean(SINR,1)));
SEmrc = sum(SEmrc,1)



SINR = zeros([N L numel(SNR)]);
Beta_l = zeros([N L]);
phi_c = zeros([1 C]);
Delta_cs = zeros([1 S]);
phi_cs = zeros([C S]);
h_l = zeros([L C S M]);
for snr = 1:numel(SNR)
    for iter = 1:N
        %determine H matrix for L users
        
          for l = 1:L  
            Mags =sqrt(r_o^2+ +(r^2-r_o^2)*rand(1,1));
            Phases = -r +(2*r)*rand(1,1);
            d_l = (Mags .* exp(1i*Phases));
            X_l = (10.^(normrnd(0,sigma_sf)/10));
            Beta_l(iter,l) = A*X_l*(abs(d_l)/d_o).^-Gamma;
            H = 1/sqrt(2)*(randn(M,L)+1i*randn(M,L));
            H = sqrt(R).*H;
          end
            
        
        
        %determine SINR for each user
        for l = 1:L
            
            W = H*inv(H'*H);
            Sum = 0;
            for lb = 1:L
                if lb ~= l
                    temp = abs(W(:,l)'*H(:,lb))^2;
                    Sum = Sum + temp;
                end
            end
            %temp = (p_ul*(abs(W(:,l)'*H(:,l)))^2)/(p_ul*Sum + (p_ul*Beta_l(iter,l)./SNRpow(snr)) *norm(W(:,l))^2);
            SINR(iter,l,snr) = (p_ul*(abs(W(:,l)'*H(:,l)))^2)/(p_ul*Sum + (p_ul*Beta_l(iter,l)./SNRpow(snr)) *norm(W(:,l))^2);
            
        end
    end
end

SEzf = squeeze(log2(1+mean(SINR,1)));
SEzf = sum(SEzf,1)

%figure(4)
%hold on
plot(SNR,(SEmrc),'r--')
plot(SNR,(SEzf),'b--')
grid on
xlabel('SNR dB')
ylabel('Spectral Efficiency bps/HZ')
title('Up-link Performace for Kronecker and Ray Based Channel')
legend('MRC Ray-Based','ZF Ray-Based','MRC Kronecker','ZF Kronecker')

%% CDFs
figure(5)
hold on
count = 1;
N = 2000;
SNR = 0;
for M = [25 100 150]
%SNR = 10;
SNRpow = db2pow(SNR)
SINR = zeros([N L numel(SNR)]);
Beta_l = zeros([N L]);
phi_c = zeros([1 C]);
Delta_cs = zeros([1 S]);
phi_cs = zeros([C S]);
h_l = zeros([L C S M]);
for snr = 1:numel(SNR)
    for iter = 1:N
        %determine H matrix for L users
        for l = 1:L
            Mags =sqrt(r_o^2+ +(r^2-r_o^2)*rand(1,1));
            Phases = -r +(2*r)*rand(1,1);
            d_l = (Mags .* exp(1i*Phases));
            X_l = (10.^(normrnd(0,sigma_sf)/10));
            Beta_l(iter,l) = A*X_l*(abs(d_l)/d_o).^-Gamma;
            phi_c = normrnd(0,sigma_c,C,1);
            
            Betac = Beta_l(iter,l)./(2.^(1:C));
            Betac(C) = Betac(C-1);
            
            for c = 1:C
                Delta_cs = laprnd(S,1,0,sigma_s);
                phi_cs(c,:) = phi_c(c) + Delta_cs;
                Beta_cs_l = Betac(c)/S;
                for s = 1:S               
                    theta_cs_l = (2*pi)*rand(1);
                    gamma_l = sqrt(Beta_cs_l)*exp(i*theta_cs_l);
                    a_l = exp(i*2*pi*(0:M-1)*delta*sind(phi_cs(c,s)))';
                    h_l(l,c,s,:) = gamma_l*a_l;                  
                end
            end
            h_l_sum = sum(h_l,2);
            h_l_sum = squeeze(sum(h_l_sum,3))';                 
        end     
        for l = 1:L            
            H = h_l_sum;
            W = H;           
            Sum = 0;
            for lb = 1:L
                if lb ~= l
                    temp = abs(W(:,l)'*H(:,lb))^2;
                    Sum = Sum + temp;
                end
            end

            SINR(iter,l,snr) = (p_ul*(abs(W(:,l)'*H(:,l)))^2)/(p_ul*Sum + (p_ul*Beta_l(iter,l)./SNRpow(snr)) *norm(W(:,l))^2);
        end
    end  
end

SEmrc = ((SINR));
SEmrc_cdf = mean(SEmrc,2);

cdfData150(:,count) = pow2db(SEmrc_cdf);

%cdfplot(SEmrc_cdf)
plots(count) = cdfplot(pow2db(SEmrc_cdf));
count = count+1;
end

for M = [25 100 150]
%SNR = 10;
SNRpow = db2pow(SNR)
SINR = zeros([N L numel(SNR)]);
Beta_l = zeros([N L]);
phi_c = zeros([1 C]);
Delta_cs = zeros([1 S]);
phi_cs = zeros([C S]);
h_l = zeros([L C S M]);
for snr = 1:numel(SNR)
    for iter = 1:N
        %determine H matrix for L users
        for l = 1:L
            Mags =sqrt(r_o^2+ +(r^2-r_o^2)*rand(1,1));
            Phases = -r +(2*r)*rand(1,1);
            d_l = (Mags .* exp(1i*Phases));
            X_l = (10.^(normrnd(0,sigma_sf)/10));
            Beta_l(iter,l) = A*X_l*(abs(d_l)/d_o).^-Gamma;
            phi_c = normrnd(0,sigma_c,C,1);
            
            Betac = Beta_l(iter,l)./(2.^(1:C));
            Betac(C) = Betac(C-1);
            
            for c = 1:C
                Delta_cs = laprnd(S,1,0,sigma_s);
                phi_cs(c,:) = phi_c(c) + Delta_cs;
                Beta_cs_l = Betac(c)/S;
                for s = 1:S               
                    theta_cs_l = (2*pi)*rand(1);
                    gamma_l = sqrt(Beta_cs_l)*exp(i*theta_cs_l);
                    a_l = exp(i*2*pi*(0:M-1)*delta*sind(phi_cs(c,s)))';
                    h_l(l,c,s,:) = gamma_l*a_l;                  
                end
            end
            h_l_sum = sum(h_l,2);
            h_l_sum = squeeze(sum(h_l_sum,3))';                 
        end     
        %determine SINR for each user
        for l = 1:L            
            H = h_l_sum;
            W = H*inv(H'*H);           
            Sum = 0;
            for lb = 1:L
                if lb ~= l
                    temp = abs(W(:,l)'*H(:,lb))^2;
                    Sum = Sum + temp;
                end
            end

            SINR(iter,l,snr) = (p_ul*(abs(W(:,l)'*H(:,l)))^2)/(p_ul*Sum + (p_ul*Beta_l(iter,l)./SNRpow(snr)) *norm(W(:,l))^2);
        end
    end  
end

SEzf = ((SINR));
SEzf_cdf = mean(SEzf,2);

cdfData150(:,count) = pow2db(SEzf_cdf);


%cdfplot(SEmrc_cdf)
plots(count) = cdfplot(pow2db(SEzf_cdf));
count = count+1;
end


set( plots(1), 'LineStyle', ':', 'Color', 'k');
set( plots(2), 'LineStyle', '--', 'Color', 'k');
set( plots(3), 'LineStyle', '-', 'Color', 'k');
set( plots(4), 'LineStyle', ':', 'Color', 'b');
set( plots(5), 'LineStyle', '--', 'Color', 'b');
set( plots(6), 'LineStyle', '-', 'Color', 'b');
xlabel('SINR dB')
ylabel('Probability')
title('SINR for Ray-Based Models at SNR = 10dB, M = 25, 100, 150')
legend('MRC M = 25','MRC M = 100','MRC M = 150','ZF M = 25','ZF M = 100','ZF M = 150')




figure (6)
hold on
for M = [25 150 500]
p = 0.9
R = zeros([M L]);
for j = 1:M
    for k = 1:L
        R(j,k) = p^abs(j-k);
    end
end 
SINR = zeros([N L numel(SNR)]);
Beta_l = zeros([N L]);
phi_c = zeros([1 C]);
Delta_cs = zeros([1 S]);
phi_cs = zeros([C S]);
h_l = zeros([L C S M]);
for snr = 1:numel(SNR)
    for iter = 1:N
        %determine H matrix for L users
        for l = 1:L
            
            Mags =sqrt(r_o^2+ +(r^2-r_o^2)*rand(1,1));
            Phases = -r +(2*r)*rand(1,1);
            d_l = (Mags .* exp(1i*Phases));
            X_l = (10.^(normrnd(0,sigma_sf)/10));
            Beta_l(iter,l) = A*X_l*(abs(d_l)/d_o).^-Gamma;
            H = 1/sqrt(2) *(randn(M,L)+1i*randn(M,L));
            H = sqrt(R).*H;
            
        end  
        
        
        %determine SINR for each user
        for l = 1:L
            
            W = H;
            Sum = 0;
            for lb = 1:L
                if lb ~= l
                    temp = abs(W(:,l)'*H(:,lb))^2;
                    Sum = Sum + temp;
                end
            end
            
            SINR(iter,l,snr) = (p_ul*(abs(W(:,l)'*H(:,l)))^2)/(p_ul*Sum + (p_ul*Beta_l(iter,l)./SNRpow(snr)) *norm(W(:,l))^2);
            
        end
    end
    
end

SEmrc = ((SINR));
SEmrc_cdf = mean(SEmrc,2);
cdfData150(:,count) = pow2db(SEmrc_cdf);

%cdfplot(SEmrc_cdf)
plots(count) = cdfplot(pow2db(SEmrc_cdf));
count = count+1;
end

for M = [25 150 500]
p = 0.9
R = zeros([M L]);
for j = 1:M
    for k = 1:L
        R(j,k) = p^abs(j-k);
    end
end 
SINR = zeros([N L numel(SNR)]);
Beta_l = zeros([N L]);
phi_c = zeros([1 C]);
Delta_cs = zeros([1 S]);
phi_cs = zeros([C S]);
h_l = zeros([L C S M]);
for snr = 1:numel(SNR)
    for iter = 1:N
        %determine H matrix for L users
        for l = 1:L
            
            Mags =sqrt(r_o^2+ +(r^2-r_o^2)*rand(1,1));
            Phases = -r +(2*r)*rand(1,1);
            d_l = (Mags .* exp(1i*Phases));
            X_l = (10.^(normrnd(0,sigma_sf)/10));
            Beta_l(iter,l) = A*X_l*(abs(d_l)/d_o).^-Gamma;
            H = 1/sqrt(2) *(randn(M,L)+1i*randn(M,L));
            H = sqrt(R).*H;
            
        end  
        
        
        %determine SINR for each user
        for l = 1:L
            
            W = H*inv(H'*H);
            Sum = 0;
            for lb = 1:L
                if lb ~= l
                    temp = abs(W(:,l)'*H(:,lb))^2;
                    Sum = Sum + temp;
                end
            end
            
            SINR(iter,l,snr) = (p_ul*(abs(W(:,l)'*H(:,l)))^2)/(p_ul*Sum + (p_ul*Beta_l(iter,l)./SNRpow(snr)) *norm(W(:,l))^2);
            
        end
    end
    p_ul*Sum
    (p_ul*Beta_l(iter,l)./SNRpow(snr))
    norm(W(:,l))^2
end

SEzf = ((SINR));
SEzf_cdf = mean(SEzf,2);
cdfData150(:,count) = pow2db(SEzf_cdf);
%cdfplot(SEmrc_cdf)
plots(count) = cdfplot(pow2db(SEzf_cdf));
count = count+1;
end

set( plots(7), 'LineStyle', ':', 'Color', 'k');
set( plots(8), 'LineStyle', '--', 'Color', 'k');
set( plots(9), 'LineStyle', '-', 'Color', 'k');
set( plots(10), 'LineStyle', ':', 'Color', 'b');
set( plots(11), 'LineStyle', '--', 'Color', 'b');
set( plots(12), 'LineStyle', '-', 'Color', 'b');
xlim([-5 15])
xlabel('SINR dB')
ylabel('Probability')
title('SINR for Kronecker Models at SNR = 10dB, M = 25, 150, 500')
legend('MRC M = 25','MRC M = 150','MRC M = 500','ZF M = 25','ZF M = 150','ZF M = 500')
%%

figure(7)
hold on
for inc = [3 6 8 11]
plots(count) = cdfplot(cdfData150(:,inc));
count =count +1;
end
set( plots(13), 'LineStyle', '--', 'Color', 'r');
set( plots(14), 'LineStyle', '-', 'Color', 'r');
set( plots(15), 'LineStyle', '--', 'Color', 'b');
set( plots(16), 'LineStyle', '-', 'Color', 'b');
xlabel('SINR dB')
ylabel('Probability')
title('SINR for all Models at SNR = 10dB, M = 150')
legend('MRC Ray-Based','ZF Ray-Based','MRC Kronecker','ZF Kronecker')
xlim([-5 25])













