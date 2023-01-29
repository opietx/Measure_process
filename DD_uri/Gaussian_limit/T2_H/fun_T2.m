T2 = zeros(length(cTemp),1);
w = zeros(length(cTemp), TimeRes);

for Temp_ind = 1:length(cTemp)
    Temp = cTemp(Temp_ind);

    Exp = alpha(Temp);
    %Exp=1.2;
    
    noiseVar = noise(Temp);
    
    for tCount = 1:TimeRes
        t = Time(tCount);
        if method == 'H'
            %%%%%%%%%%%%%%%%%%%          HAHN ECHO        %%%%%%%%%%%%%%%%%%%%%%%
            fun = @(AngFreq) (2*noiseVar/(pi*wNorm))*((atan(2*rMax./AngFreq)- atan(2*rMin./AngFreq))./AngFreq.^Exp).* ...
                        (8*(sin(AngFreq*t/4).^4))./ ...
                        (AngFreq.^2);        
        elseif method == 'S'
            %%%%%%%%%%%%%%%%%%%          FID        %%%%%%%%%%%%%%%%%%%%%%%
            fun = @(AngFreq) (2*noiseVar/(pi*wNorm))*((atan(2*rMax./AngFreq)- atan(2*rMin./AngFreq))./AngFreq.^Exp).* ...
                        (2*(sin(AngFreq*t/2).^2))./ ...
                        (AngFreq.^2);     
        end
        chi = integral(fun, 1e-5, 1e10);
        w(Temp_ind,tCount) = exp(-chi);

    end        
    [~, IndexT2] = min(abs(w(Temp_ind,:) - exp(-1)));
    T2(Temp_ind) = Time(IndexT2);

    
end
