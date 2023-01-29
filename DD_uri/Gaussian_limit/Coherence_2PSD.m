clear
clc
tic

%# import all params used in the experiment
params;

method = 'CPMG';

if method == 'S'
    test
    plot_T2
elseif method == 'H'
    test
    plot_T2
%     f = fit(cTemp',T2*1e6,'power2');
%     plot(f)
%     params_fit = coeffvalues(f);
%     params_fit(2)
elseif method == 'CPMG'
    cNPi = [2 4 8 16 20];
    test2
    plot_DD
end

toc



