% Temperature-dependence of the Ramsey decay for a set of qubits embeded in
% inhomogeneous baths of tunneling Two-Level Systems (TLSs). Following the
% common convention we set hbar = 1. 

tic        % Cornometer start
clc       % Clear workspace editor
clear    % Erase workspace memory
 
% %%%     Model Constants     %%%
nrTLSs = 1000;       % Number of Two-Level Systems (TLSs)
nrQubits = 1;      % Number of Qubits
engMin = 0.005;   % Minimum TLS energy. This is the energy difference between the minima of the 
                             % double-well potential. This energy can be used
                             % as a measure of TLS asymettry. Units degree Kelvin.
engMax = 0.02;     % Maximum TLS energy. Units degree Kelvin.
resTemp = 100;     % Resolution of temperature.
minTemp = .01;     % Minimum temperature. Units ueV.
maxTemp = 2;     % Maximum temperature. Units ueV.
 

asymEngMat = engMin + (engMax - engMin)*rand(nrQubits, nrTLSs);
engMin2 = 0.01;   % Minimum TLS energy. This is the energy difference between the minima of the 
engMax2 = 0.01;      % Maximum TLS energy. Units degree Kelvin.

asymEngMat2 = engMin2 + (engMax2 - engMin2)*rand(nrQubits, nrTLSs);

vCoupDummy = 0.1 + 0.01*randn(nrQubits, nrTLSs);
vCoupMat = exp(vCoupDummy*1e6);               % Normally distributed coupling strengths (mean 4e4 s^{-1} and deviation 1e4 s^{-1}).
vCoupDummy2 = 5 + .3*randn(nrQubits, nrTLSs);
vCoupMat2 = exp(vCoupDummy2);               % Normally distributed coupling strengths (mean 4e4 s^{-1} and deviation 1e4 s^{-1}).


temp = linspace(minTemp, maxTemp, resTemp);                                   % List of temperatures
freqMod = zeros(nrQubits, resTemp);



for qubitCount = 1:nrQubits
    freqModTot = zeros(1, resTemp);
   
    for TLSCount = 1:nrTLSs
        vCoup = vCoupMat(qubitCount, TLSCount);
        vCoup2 = vCoupMat2(qubitCount, TLSCount);
                
        asymEng = asymEngMat(qubitCount, TLSCount);
        asymEng2 = asymEngMat2(qubitCount, TLSCount);
        
        
        rateExcite = 2*vCoup.*exp(-asymEng./(2*temp))./sinh(asymEng./(2*temp));
        rateRelax = 2*vCoup.*exp(asymEng./(2*temp))./sinh(asymEng./(2*temp));
       
        rate_m = rateExcite-rateRelax;
        rate_p = rateExcite+rateRelax;
        freqModSingle = rate_m./rate_p;
%         freqModSingle = -vCoup*tanh(asymEng./temp);
        freqModTot = freqModTot + freqModSingle;
    end

    freqMod(qubitCount, :) = freqModTot - freqModTot(1, 1);

end


box on

set(gca,'LineWidth',1,'TickLabelInterpreter','latex','FontSize',24)
xlabel('$T~(^\circ \rm{K})$','interpreter','latex','FontSize',24)
ylabel('$\delta E_Z~(\rm{Hz})$','interpreter','latex','FontSize',24)

for qubitCnt = 1:nrQubits
    
    hold on
    plot(temp, freqMod(qubitCnt, :), 'LineWidth', 3)     
    
end

toc    % Cornometer stop
