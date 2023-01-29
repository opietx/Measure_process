clear
clc
tic

TimeRes = 300;
Time = logspace(-8, -3, TimeRes);
rMin = 1e-3;
rMax = 1e8;
wNorm = log(rMax/rMin);


method = 'CPMG';

base_noise = 1.5e11;
cTemp = [200,300,400,500,600,700,800,900,1000];
alpha = @(temperature) 1.25-0.95*temperature./1000;
noise = @(temperature) base_noise*(1 + 0.5*(temperature./1000)^2);
noiseVar = base_noise;
if method == 'S'
    fun_T2
    plot_T2
elseif method == 'H'
    fun_T2
    plot_T2
elseif method == 'CPMG'
    cNPi = [2 8 20];
    fun_DD
    plot_DD
end

toc
