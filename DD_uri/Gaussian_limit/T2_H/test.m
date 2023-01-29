T2 = zeros(length(cTemp),1);
w = zeros(length(cTemp), TimeRes);

for Temp_ind = 1:length(cTemp)
    Temp = cTemp(Temp_ind)

    noise1 = noise_f(Temp);
    noise2 = noise_w(Temp);
    
    for tCount = 1:TimeRes
        t = Time(tCount);
        if method == 'H'
            %%%%%%%%%%%%%%%%%%%          HAHN ECHO        %%%%%%%%%%%%%%%%%%%%%%%
            fun = @(AngFreq) (2*noiseVar/(pi*wNorm))*(atan(2*rMax./AngFreq)- atan(2*rMin./AngFreq)).*(noise1./(AngFreq.^alpha1) + noise2./(AngFreq.^alpha2)).* ...
                        (8*(sin(AngFreq*t/4).^4))./ ...
                        (AngFreq.^2);        
        elseif method == 'S'
            %%%%%%%%%%%%%%%%%%%          FID        %%%%%%%%%%%%%%%%%%%%%%%
            fun = @(AngFreq) (2*noiseVar./(pi*wNorm))*(atan(2*rMax./AngFreq)- atan(2*rMin./AngFreq)).*(noise1./AngFreq.^alpha1 + noise2./AngFreq.^alpha2).* ...
                        (2*(sin(AngFreq*t/2).^2))./ ...
                        (AngFreq.^2);     
        end
        chi = integral(fun, 1e-8, 1e10);
        w(Temp_ind,tCount) = exp(-chi);

    end        
    [~, IndexT2] = min(abs(w(Temp_ind,:) - exp(-1)));
    T2(Temp_ind) = Time(IndexT2);

    
end
