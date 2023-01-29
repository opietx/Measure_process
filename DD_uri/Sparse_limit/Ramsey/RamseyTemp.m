
  1
  2
  3
  4
  5
  6
  7
  8
  9
 10
 11
 12
 13
 14
 15
 16
 17
 18
 19
 20
 21
 22
 23
 24
 25
 26
 27
 28
 29
 30
 31
 32
 33
 34
 35
 36
 37
 38
 39
 40
 41
 42
 43
 44
 45
 46
 47
 48
 49
 50
 51
 52
 53
 54
 55
 56
 57
 58
 59
 60
 61
 62
 63
 64
 65
 66
 67
 68
 69
 70
 71
 72
 73
 74
 75
 76
 77
 78
 79
 80
 81
 82
 83
 84
 85
 86
 87
 88
 89
 90
 91
 92
 93
 94
 95
 96
 97
 98
 99
100
101
102
103
104
105
106
107
108
% Temperature-dependence of the Ramsey decay for a set of qubits embeded in
% inhomogeneous baths of tunneling Two-Level Systems (TLSs). Following the
% common convention we set hbar = 1. 

tic        % Cornometer start
clc       % Clear workspace editor
clear    % Erase workspace memory
 
% %%%     Model Constants     %%%
nrTLSs = 1e3;       % Number of Two-Level Systems (TLSs)
nrQubits = 3;      % Number of Qubits
rMin = 1e-3;            % Minimum TLS jumping rate. Units s^{-1}.
rMax = 1e6;        % Maximum TLS jumping rate. Units s^{-1}
engMin = .01;   % Minimum TLS energy. This is the energy difference between the minima of the 
                             % double-well potential. This energy can be used
                             % as a measure of TLS asymettry. Units degree Kelvin.
engMax = 100;      % Maximum TLS energy. Units degree Kelvin.
finTime = 1e-2;   % Final time for simulated Ramsey decay. Units seconds.
resTime = 1e4;    % Resolution of time.
resTemp = 10;     % Resolution of temperature.
minTemp = .1;     % Minimum temperature. Units ueV.
maxTemp = 1.5;     % Maximum temperature. Units ueV.

% %%%     Variab1.5e ranges     %%%
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

for qubitCount = 1:nrQubits
    
    for tempCount = 1:resTemp
        decTot = ones(1, length(time));
        tempBoltz = tempList(tempCount)
        rateExciteMat = 2*rateBare.*exp(-asymEngMat/(2*tempBoltz))./sinh(asymEngMat/(2*tempBoltz));
        rateRelaxMat = 2*rateBare.*exp(asymEngMat/(2*tempBoltz))./sinh(asymEngMat/(2*tempBoltz));
        
            for TLSCount = 1:nrTLSs

                rateExcite = rateExciteMat(qubitCount, TLSCount);
                rateRelax = rateRelaxMat(qubitCount, TLSCount);
                vCoup = vCoupMat(qubitCount, TLSCount);
                rateDiff = rateRelax - rateExcite;
                rateSum = rateExcite + rateRelax;
                Gamma = -1i*vCoup*rateDiff/rateSum + rateSum;
                Alpha = sqrt(-vCoup^2 - 2i*vCoup*rateDiff + rateSum^2);

                RamseySingle =.5*real((1+Gamma/Alpha)*exp((-Gamma + Alpha)*time) + (1-Gamma/Alpha)*exp((-Gamma - Alpha)*time)) ;
                
                decTot = decTot.*RamseySingle;
                
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

 for qubitCount = 1:nrQubits
     for tempCount = 1:resTemp
         hold on
         colorRGB = ColorMat(qubitCount, :);
         plot(tempList(tempCount), T2starMat(qubitCount, tempCount), ...
             'o', 'MarkerSize', 10, 'Color', 'k', 'LineWidth', 1, 'MarkerFaceColor', ColorMat(qubitCount,:))
     end
 end

toc    % Cornometer stop
