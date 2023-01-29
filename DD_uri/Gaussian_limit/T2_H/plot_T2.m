%%%%%%%%%%%%%%FIGURE WITH DECAYS%%%%%%%%%%%%%%%%%%%%%%%%%
fig1 =figure(1);
box on
hold on
    
for Temp_ind = 1:length(cTemp)
    plot(Time, w(Temp_ind,:), '-', ...
        'LineWidth',2,'Color','k','Color',[Temp_ind/length(cTemp) .1 .1])
    
end
set(gca, 'LineWidth', 2)
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 16)
xlabel('Time ~[arb. units]','Interpreter','LaTeX','FontSize',22)
ylabel('Spin prob. ~[arb. units]','Interpreter','LaTeX','FontSize',22)
title(method);
set(gca,'xscale','log')





%%%%%%%%%%%%%%FIGURE WITH TIME EVOLUTION%%%%%%%%%%%%%%%%%%%%%%%%%
fig2 =figure(2);
box on
hold on
    
for Temp_ind = 1:length(cTemp)
    plot(cTemp(Temp_ind), T2(Temp_ind)*1e6, 'o', ...
        'LineWidth',2,'MarkerSize',14,'Color','k','MarkerFaceColor',[Temp_ind/length(cTemp) .1 .0])
    
end

set(gca, 'LineWidth', 2)
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 16)


xlabel('Temperature [mK]','Interpreter','LaTeX','FontSize',22)
ylabel('$T_2^{\rm{H}}$~[$\mu$s]','Interpreter','LaTeX','FontSize',22)



title(method);

toc
