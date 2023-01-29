clc
clear

tempBoltz = 16;

nrTLSs = 1000;
rMax = 10^3;
rMin = 10^-3;
freq = logspace(log10(rMin) -10, log10(rMax) + 10, 2000);


engMin = 0.0001;
engMax = 0.01;

vCoupling = 1;


specTot = zeros(1,length(freq));

param = tempBoltz;
name = 'temp';


factor = 1;
for tlsCount = 1:nrTLSs
       rateBare =  rMin*exp(log(rMax/rMin)*rand);
       e1 = engMin + (engMax-engMin)*rand;      
       rateExcite = 2*rateBare.*exp(-e1/(2*tempBoltz));
       
       
       e2 = engMin + (engMax-engMin)*rand;      
       rateRelax = 2*rateBare.*exp(-e2/(2*tempBoltz));
       
       
       rate = rateExcite + rateRelax;
    
       vCoup = vCoupling*rand;
       
       specSingle = 4*(vCoup^2)*(1/rateRelax+1/rateExcite)./(4*(rateRelax+rateExcite)^2 + freq.^2);
       specTot = specTot + specSingle;
       
    
end
plot_PSD