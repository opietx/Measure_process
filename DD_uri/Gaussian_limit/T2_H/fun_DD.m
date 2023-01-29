T2 = zeros(length(cTemp),length(cNPi));
w = zeros(length(cTemp), TimeRes, length(cNPi));


for Temp_ind = 1:length(cTemp)

    Temp = cTemp(Temp_ind);

    Exp = alpha(Temp)
    %Exp=1.2;
    %noiseVar = noise(Temp);

    for nCyc = 1:length(cNPi)
    
        nPi = cNPi(nCyc)

        for tCount = 1:TimeRes   
            t = Time(tCount);
         
            fun = @(AngFreq) (2*noiseVar/(pi*wNorm))*((atan(2*rMax./AngFreq) - atan(2*rMin./AngFreq))./AngFreq.^Exp).*...
                (8*(sin(AngFreq*t/(4*nPi)).^4).*(sin(AngFreq*t/2).^2)./(cos(AngFreq*t/(2*nPi)).^2))./ ...
                (AngFreq.^2);
            
            chi = integral(fun, 1e-6, 1e10);
            w(Temp_ind,tCount,nCyc) = exp(-chi);
  
        end
        [~, IndexT2] = min(abs(w(Temp_ind,:,nCyc) - exp(-1)));
        T2(Temp_ind,nCyc) = Time(IndexT2);
    end
end