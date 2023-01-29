%# import all params used in the experiment
% clear
% clc
tic
% params;


method = 'S';
test;
T2s = T2;

method = 'H';
test;
T2h = T2;


method = 'C';
T2 = T2h./T2s;
fig1 = figure;
box on
hold on
    
for Temp_ind = 1:length(cTemp)
    plot(cTemp(Temp_ind), T2(Temp_ind), 'o', ...
        'LineWidth',2,'MarkerSize',14,'Color','k','MarkerFaceColor',[Temp_ind/length(cTemp) .1 .0])
    
end

set(gca, 'LineWidth', 2)
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 16)


xlabel('Temperature [mK]','Interpreter','LaTeX','FontSize',22)
ylabel('$T_2^{\rm{H}}/T_2^{\rm{*}}$','Interpreter','LaTeX','FontSize',22)


toc