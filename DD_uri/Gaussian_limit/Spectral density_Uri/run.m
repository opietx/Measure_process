clc
clear

tempBoltz = 150;

nrTLSs = 1;
rMax = 10^8;
rMin = 10^-3;
freq = logspace(log10(rMin) -10, log10(rMax) + 5, 2000);


engMin = 80;
engMax = 5;

vCoupling = 1;


specTot = zeros(1,length(freq));

param = tempBoltz;
name = 'temp';


factor = 1
for tlsCount = 1:nrTLSs
       rateBare =          0.00001;%rMin*exp(log(rMax/rMin)*rand);
       
       e1 = abs(normrnd(engMin,engMax));
       e2 = abs(normrnd(engMin*factor,engMax*factor));
       asymEngMat = 8000;%(e1^2+e2^2)^0.5;
       
       rateExcite = 2*rateBare.*exp(-asymEngMat/(2*tempBoltz));%./sinh(asymEngMat/(2*tempBoltz));
       rateRelax = 2*rateBare.*exp(asymEngMat/(2*tempBoltz));%./sinh(asymEngMat/(2*tempBoltz));
       rate = rateExcite + rateRelax;

%        rate = rMin*exp(log(rMax/rMin)*rand);
    
       vCoup = vCoupling;%*rand;
       
       specSingle = 4*(vCoup^2)*rate./(4*rate^2 + freq.^2);
       specTot = specTot + specSingle;
       
    
end
plot_PSD