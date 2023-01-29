
figure
hold on
Omega = logspace(-3,8,200);
rMin = 1e-5;
rMax = 1e9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%COLD%%%%%%%%%%%%%%%%%%%%%%%%%%

window_1_f = 0.01*(atan(2*rMax./(Omega)) - atan(2*rMin./(Omega)))./Omega.^1.2;
window = window_1_f;

plot(Omega, window, '-','LineWidth', 3, 'Color', [.1 .0 .9], 'DisplayName','Cold');




set(gca, 'LineWidth', 2);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 16);
xlabel('$\omega~(\rm{s}^{-1})$', 'Interpreter', 'LaTeX', 'FontSize', 22);
ylabel('PSD [arb. units]','Interpreter','LaTeX','FontSize',22);
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

lims = ylim
patch([1e-1 20 20 1e-1], [max(lims)*[1 1] min(lims)*[1,1] ], [0.8 0.8 0.8]);
patch([1e2 2e2 2e2 1e2], [max(lims)*[1 1] min(lims)*[1,1] ], [0.8 0.8 0.8]);
patch([14e3 1e5 1e5 14e3], [max(lims)*[1 1] min(lims)*[1,1] ], [0.8 0.8 0.8]);
set(gca,'children',flipud(get(gca,'children')));