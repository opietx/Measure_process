clear
clc

Omega = logspace(-2,7,600)*(2*pi);
params;

factor = 3e-3;
save =0;

figure
hold on
for temp_ind = 1:length(cTemp)
    Temp = cTemp(temp_ind);
    noise1 = noise_f(Temp);
    noise2 = noise_w(Temp);
    
    

    window_1_f = noise1*(2*noiseVar/(pi*wNorm))*(atan(2*rMax./(Omega)) - atan(2*rMin./(Omega)))./Omega.^alpha1;
    window_white = noise2*(2*noiseVar/(pi*wNorm))*(atan(2*rMax./(Omega)) - atan(2*rMin./(Omega)))./Omega.^alpha2;

    window = (window_1_f+window_white)*factor;


%     plot(Omega, window_1_f,'--', 'LineWidth', 1, 'Color', [.1 .0 .1])
%     plot(Omega, window_white,'--', 'LineWidth', 1, 'Color', [.1 .0 .1])
    plot(Omega, window, '-','LineWidth', 3, 'Color', [temp_ind/length(cTemp) .0 .1], 'DisplayName',string(Temp));
    if save==1
        m = [Omega', window'];
        writematrix(m,'M:\tnw\ns\qt\spin-qubits\projects\Hot qubit SiGe\Data\T2_CPMG_simulation\data\PSD\'+string(Temp)+'mK.csv')
    end
    
    
end


set(gca, 'LineWidth', 2);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 16);
xlabel('$\omega~(\rm{s}^{-1})$', 'Interpreter', 'LaTeX', 'FontSize', 22);
ylabel('PSD [arb. units]','Interpreter','LaTeX','FontSize',22);
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

patch([1e-1 20 20 1e-1], [max(ylim)*[1 1] min(ylim)*[1,1] ], [0.8 0.8 0.8]); %SET
patch([14.8e3 256e3 256e3 14.8e3], [max(ylim)*[1 1] min(ylim)*[1,1] ], [0.8 0.8 0.8]);
% legend;
xlim([1e-2, 1e7])
ylim([1-2, 1e10])

set(gca,'children',flipud(get(gca,'children')));


