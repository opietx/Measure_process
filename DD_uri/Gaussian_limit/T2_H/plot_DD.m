%%%%%%%%%%%%%%FIGURE WITH DECAYS%%%%%%%%%%%%%%%%%%%%%%%%%
fig1 = figure(1);
box on
hold on
    
for Temp_ind = 1:length(cTemp)
    for nCyc = 1:length(cNPi)
        plot(Time, w(Temp_ind,:, nCyc), '-', ...
            'LineWidth',2,'Color','k','Color',[Temp_ind/length(cTemp) .1 nCyc/length(cNPi)])
    end
end

set(gca, 'LineWidth', 2)
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 16)
b_bb = gca; legend(b_bb,'off');
xlabel('Time ~[arb. units]','Interpreter','LaTeX','FontSize',22)
ylabel('Spin prob. ~[arb. units]','Interpreter','LaTeX','FontSize',22)
title(method);
set(gca,'xscale','log')



hold off

%%%%%%%%%%%%%%FIGURE WITH TIME EVOLUTION%%%%%%%%%%%%%%%%%%%%%%%%%
fig2 = figure(2);
box on
hold on
    
for Temp_ind = 1:length(cTemp)
    for nCyc = 1:length(cNPi)    
        plot(cTemp(Temp_ind), T2(Temp_ind,nCyc)*1e6, 'o', ...
            'LineWidth',2,'MarkerSize',14,'Color','k','MarkerFaceColor',[Temp_ind/length(cTemp) .1 nCyc/length(cNPi)])
    end
end

set(gca, 'LineWidth', 2)
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 16)
set(gca, 'xscale', 'linear')
set(gca, 'yscale', 'log')
b_bb = gca; legend(b_bb,'off');




xlabel('Temperature [mK]','Interpreter','LaTeX','FontSize',22)
ylabel('$T_2^{\rm{CPMG}}$~[$\mu$s]','Interpreter','LaTeX','FontSize',22)
title(method);

hold off
%%%%%%%%%%%%%%%%%%% figure with N-pi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig3 = figure(3);
box on
hold on

alphas = zeros(1,length(cTemp));
for Temp_ind = 1:length(cTemp)
     
        plot(cNPi, T2(Temp_ind,:)*1e6, 'o', ...
            'LineWidth',2,'MarkerSize',14,'Color','k','MarkerFaceColor',[Temp_ind/length(cTemp) .1 0.1])
        f = fit(cNPi',T2(Temp_ind,:)'*1e6,'power1');
        plot(f, 'k--');
        params_fit = coeffvalues(f);
        alphas(Temp_ind) = params_fit(2)./(1-params_fit(2));
        
        
    
end


set(gca, 'LineWidth', 2)
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 16)
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')

xlabel('Number of $\pi$ pulses $n_\pi$','Interpreter','LaTeX','FontSize',22)
ylabel('$T_2^{\rm{CPMG}}$~[$\mu$s]','Interpreter','LaTeX','FontSize',22)
xlim([1,30])
b_bb = gca; legend(b_bb,'off');
hold off


%%%%%%%%%%%%%%%%%%% figure with N-pi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig4 = figure(4);
box on
hold on
plot(cTemp, alphas, 'o-', ...
            'LineWidth',2,'MarkerSize',14,'Color','k','MarkerFaceColor',[Temp_ind/length(cTemp) .1 0.1])
        

set(gca, 'LineWidth', 2)
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 16)
b_bb = gca; legend(b_bb,'off');

xlabel('Temperature [mK]','Interpreter','LaTeX','FontSize',22)
ylabel('$\alpha$','Interpreter','LaTeX','FontSize',22)
ylim([0.0,1.5])
hold off

