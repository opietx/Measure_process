clc
clear

tempBoltz =20;

nrTLSs =1e4*tempBoltz;
rMax = 10^8;
rMin = 10^-15;
freq = logspace(-2, 8 , 2000);


engMax = 100;
engMin = 0.01;

vCoupling = tempBoltz^0.5;


specTot = zeros(1,length(freq));

param = tempBoltz;
name = 'temp';

EngDis = engMin + (engMax-engMin)*rand(1,nrTLSs);%
% h = histogram(EangDis);
RateDist = rMin*exp(log(rMax/rMin)*rand(1,nrTLSs));
% h = histogram(RateDist);
factor = 1
for tlsCount = 1:nrTLSs
       rateBare = RateDist(tlsCount);
       
       asymEngMat = EngDis(tlsCount);%
       
       rateExcite = 2*rateBare.*exp(-asymEngMat/(2*tempBoltz))./sinh(asymEngMat/(2*tempBoltz));
       rateRelax = 2*rateBare.*exp(asymEngMat/(2*tempBoltz))./sinh(asymEngMat/(2*tempBoltz));
       rate = rateExcite + rateRelax;
       
%        rate = rMin*exp(log(rMax/rMin)*rand);
    
       vCoup = vCoupling;%*rand;
       if rate > rMax/100
           rate = rMax/100;
           specSingle = 4*(vCoup^2)*rate./(4*rate^2);
       else
           specSingle = 4*(vCoup^2)*rate./(4*rate^2 + freq.^2);
       end
       
       
       specTot = specTot + specSingle;
       
    
end
plot_PSD