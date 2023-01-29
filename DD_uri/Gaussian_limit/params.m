%#   time range 
TimeRes = 500;
Time = logspace(-9, -2, TimeRes);
rMin = 1e-5;
rMax = 1e10;
wNorm = log(rMax/rMin);


%# 1-f noise properties
alpha1=1.3;
alpha2=0.0;
A = 8e1; %8e1 *4
B = 80e5; %80e5
noise_f = @(temperature) (A*temperature^2) + B;
noiseVar = 0.65e5; %0.7e5

%# white-noise properties
n2_gain = 6e-6;%6e-6  
n2_base = 0.015; %0.015
noise_w = @(temperature) n2_base + n2_gain*(temperature)^2;


cTemp = linspace(200,900,5);
