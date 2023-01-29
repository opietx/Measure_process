figure
hold on

Omega = logspace(-3,8,200);
rMin = 1e-5;
rMax = 1e9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%COLD%%%%%%%%%%%%%%%%%%%%%%%%%%

window_1_f = 0.01*(atan(2*rMax./(Omega)) - atan(2*rMin./(Omega)))./Omega.^1.2;
window_white = 1e-9*(atan(2*rMax./(Omega)) - atan(2*rMin./(Omega)))./Omega.^0.1;

window = window_1_f+window_white;

plot(Omega, window_1_f,'--', 'LineWidth', 1, 'Color', [.1 .0 .3])
plot(Omega, window_white,'--', 'LineWidth', 1, 'Color', [.1 .0 .6])
cold =plot(Omega, window, '-','LineWidth', 3, 'Color', [.1 .0 .9], 'DisplayName','Cold');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%HOT%%%%%%%%%%%%%%%%%%%%%%%%%%
window_1_f = (atan(2*rMax./(Omega)) - atan(2*rMin./(Omega)))./Omega.^1.2;
window_white = 1e-5*(atan(2*rMax./(Omega)) - atan(2*rMin./(Omega)))./Omega.^0.1;

window = window_1_f+window_white;

plot(Omega, window_1_f,'--', 'LineWidth', 1, 'Color', [.3 .0 .1])
plot(Omega, window_white,'--', 'LineWidth', 1, 'Color', [.6 .0 .1])
hot = plot(Omega, window,'-', 'LineWidth', 3, 'Color', [.9 .0 .1], 'DisplayName','Hot');


set(gca, 'LineWidth', 2);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 16);
xlabel('$\omega~(\rm{s}^{-1})$', 'Interpreter', 'LaTeX', 'FontSize', 22);
ylabel('$\arctan(2\gamma_M/\omega) - \arctan(2\gamma_m/\omega)$','Interpreter','LaTeX','FontSize',22);
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

patch([1e-1 20 20 1e-1], [max(ylim)*[1 1] min(ylim)*[1,1] ], [0.8 0.8 0.8]);
patch([1e2 2e2 2e2 1e2], [max(ylim)*[1 1] min(ylim)*[1,1] ], [0.8 0.8 0.8]);
patch([14e3 1e5 1e5 14e3], [max(ylim)*[1 1] min(ylim)*[1,1] ], [0.8 0.8 0.8]);
set(gca,'children',flipud(get(gca,'children')));
l=legend([hot,cold]);


xlim([1e-4,1e9])
ylim([1e-12,1e5])
