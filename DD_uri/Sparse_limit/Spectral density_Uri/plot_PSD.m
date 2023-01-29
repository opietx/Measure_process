figure(1);
hold on

plot(freq, specTot,'LineWidth',2)
set(gca, 'LineWidth', 2)
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 16)
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')


xlabel('$\omega$ [Hz]','Interpreter','LaTeX','FontSize',22)
ylabel('PSD [arb units]','Interpreter','LaTeX','FontSize',22)
xlim([0.5e-2,1e8]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ d, ix ] = min( abs( freq-1 ) );
% val = (specTot(ix))^0.5;
% 
% figure(2);
% hold on
% 
% plot(param, val,'o', 'LineWidth',2,'Color',[0.62, 0.17, 0.41])
% set(gca, 'LineWidth', 2)
% set(gca, 'TickDir', 'out')
% set(gca, 'FontSize', 16)
% 
% title(name);
% xlabel('param [arb units]','Interpreter','LaTeX','FontSize',22)
% ylabel('PSD @ 1Hz [arb units]','Interpreter','LaTeX','FontSize',22)