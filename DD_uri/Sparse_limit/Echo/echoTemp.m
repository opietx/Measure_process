% Temperature-dependence of the Ramsey decay for a set of qubits embeded in
% inhomogeneous baths of tunneling Two-Level Systems (TLSs). Following the
% common convention we set hbar = 1. 

tic        % Cornometer start
clc       % Clear workspace editor
clear    % Erase workspace memory
 
% %%%     Model Constants     %%%
nrTLSs = 1e4;      % Number of Two-Level Systems (TLSs)
nrQubits = 1;       % Number of Qubits
rMin = 1e-1;          % Minimum TLS jumping rate. Units s^{-1}.
rMax = 1e8;         % Maximum TLS jumping rate. Units s^{-1}
engMin = .01;        % Minimum TLS energy. This is the energy difference between the minima of the 
                            % double-well potential. This energy can be used
                            % as a measure of TLS asymettry. Units degree Kelvin.
engMax = 100;    % Maximum TLS energy. Units degree Kelvin.
finTime = 1e-5;   % Final time for simulated Ramsey decay. Units seconds.
resTime = 1e2;    % Resolution of time.
resTemp = 10;     % Resolution of temperature.
minTemp = .1;     % Minimum temperature. Units ueV.
maxTemp = 1;     % Maximum temperature. Units ueV.

% %%%     Variable ranges     %%%
rateBare = rMin*exp(log(rMax/rMin)*rand(nrQubits, nrTLSs));     % TLS rate converges to a constant 
                                                                                             % asymptotically as temperature goes to infinity. 
                                                                                             % This rate is assumed to be constant.
                                                                                             % This rate is assumed to be log-uniformly distributed.
                                                                                             % [rows = qubits, columns = TLS rates]
% asymEngMat = engMin + (engMax - engMin)*rand(nrQubits, nrTLSs);   % TLS energies. [rows = qubits, columns = TLS energies]. 
%engDummy = 2.2*randn(nrQubits, nrTLSs);
% asymEngMat = exp(engDummy);
asymEngMat = engMin + (engMax - engMin)*rand(nrQubits, nrTLSs);
vCoupMat = 1e3 + .001*randn(nrQubits, nrTLSs);
tempList = linspace(minTemp, maxTemp, resTemp);                                   % List of temperatures
time = linspace(0, finTime, resTime);                                 % Time vector. 
T2starMat = zeros(nrQubits, resTemp);                  % T2stars  [rows = qubits, columns = temperatures]

ColorMat = [.15, .27, .33; 
     .16, .62, .56;
     .96, .87, .65;
     .94, .68, .39;
     1, .5, .25;
     .83, .24, .1];
 
 figure
 box on
 set(gca,'LineWidth',1,'TickLabelInterpreter','latex','FontSize',20)

for qubitCount = 1:nrQubits
    
    for tempCount = 1:resTemp
        decTot = ones(1, length(time));
        tempBoltz = tempList(tempCount)
        % rateExciteMat = 2*rateBare.*exp(-asymEngMat/(2*tempBoltz))./sinh(asymEngMat/(2*tempBoltz));
        rateExciteMat = 2*rateBare.*exp(-asymEngMat/(2*tempBoltz))./(exp(asymEngMat/(2*tempBoltz)) + exp(-asymEngMat/(2*tempBoltz))); 
        % rateRelaxMat = 2*rateBare.*exp(asymEngMat/(2*tempBoltz))./sinh(asymEngMat/(2*tempBoltz));
        rateRelaxMat = 2*rateBare.*exp(asymEngMat/(2*tempBoltz))./(exp(asymEngMat/(2*tempBoltz)) + exp(-asymEngMat/(2*tempBoltz)));
        
            for TLSCount = 1:nrTLSs

                rateExcite = rateExciteMat(qubitCount, TLSCount);
                rateRelax = rateRelaxMat(qubitCount, TLSCount);
                vCoup = vCoupMat(qubitCount, TLSCount);
                rateDiff = rateRelax - rateExcite;
                rateSum = rateExcite + rateRelax;
                Gamma = -1i*vCoup*rateDiff/rateSum + rateSum;
                Alpha = sqrt(-vCoup^2 - 2i*vCoup*rateDiff + rateSum^2);
                
                echoSingle = real(.25*exp(-rateSum*time).*((1 + (vCoup^2 + rateSum^2)/(Alpha*Alpha') + Gamma'/Alpha + Gamma/Alpha')*exp((Alpha + Alpha')*time/2) + ...
                    (1 - (vCoup^2 + rateSum^2)/(Alpha*Alpha') + Gamma'/Alpha - Gamma/Alpha')*exp((Alpha - Alpha')*time/2) + ...
                    (1 + (vCoup^2 + rateSum^2)/(Alpha*Alpha') - Gamma'/Alpha - Gamma/Alpha')*exp(-(Alpha + Alpha')*time/2) + ...
                    (1 - (vCoup^2 + rateSum^2)/(Alpha*Alpha') - Gamma'/Alpha + Gamma/Alpha')*exp((Alpha' - Alpha)*time/2)));
                
                decTot = decTot.*echoSingle;
                
            end
            
            hold on
            
            plot(time, decTot, 'LineWidth', 2, 'Color', ColorMat(qubitCount,:))
            [~, IndexT2star] = min(abs(decTot - exp(-1)));         % T2star definition and extraction.
            T2starTemp = time(IndexT2star);
            T2starMat(qubitCount, tempCount) = T2starTemp;
            
    end
    
 end

% figure('Position',[0 0 500 400])
% box on
% ylim([4.5e-6, 9.5e-6])
% yticks([5e-6 6e-6 7e-6 8e-6 9e-6])
% yticklabels({'5','6','7','8','9'})
% 
% xlim([5, 105])
% xticks([20 40 60 80 100])
% xticklabels({'20','40','60','80','100'})
% 
% set(gca,'LineWidth',1,'TickLabelInterpreter','latex','FontSize',20)
% xlabel('$k_BT~(\mu eV)$','interpreter','latex','FontSize',20)
% ylabel('$T_2^*~(\mu s)$','interpreter','latex','FontSize',20)

figure
box on
set(gca,'LineWidth',1,'TickLabelInterpreter','latex','FontSize',20)

 for qubitCount = 1:nrQubits
     for tempCount = 1:resTemp
         hold on
         colorRGB = ColorMat(qubitCount, :);
         plot(tempList(tempCount), T2starMat(qubitCount, tempCount), ...
             'o', 'MarkerSize', 15, 'Color', 'k', 'LineWidth', 1, 'MarkerFaceColor', ColorMat(qubitCount,:))
     end
 end

toc    % Cornometer stop
