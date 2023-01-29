for ii = 1:50
    
%     clear
%     clc
    params;
    ii;
    text =string(ii)+' : '+string(A)+', '+string(n2_gain)+', '+string(n2_base)+'\n';
    fid = fopen('figures/logs.txt','a');
    fprintf(fid, text);
    fclose(fid);
    
    
    method = 'CPMG';
    Coherence_2PSD
    saveas(fig4,'figures/plots/'+string(ii)+'_fig4_CPMG_.png')
    saveas(fig3,'figures/plots/'+string(ii)+'_fig3_CPMG_.png')
    
    for kk = 1:length(cTemp)
        m = [[2 4 8 16 20]', T2(kk,:)'];
        writematrix(m,'M:\tnw\ns\qt\spin-qubits\projects\Hot qubit SiGe\Data\T2_CPMG_simulation\data\CPMG\'+string(cTemp(kk))+'_CPMG.csv')
    end
    
    close all hidden
%     clear
%     clc
    
    
    method = 'H';
    Coherence_2PSD
    saveas(fig2,'figures/plots/'+string(ii)+'_fig2_H_.png')
    m = [cTemp', T2];
    writematrix(m,'figures/data/'+string(ii)+'_H_.csv')
    close all hidden
%     clear
%     clc
    
    
    method = 'S';
    Coherence_2PSD
    saveas(fig2,'figures/plots/'+string(ii)+'_fig2_S_.png')
    m = [cTemp', T2];
    writematrix(m,'figures/data/'+string(ii)+'_S_.csv')
    
    
    close all hidden
%     clear
%     clc
    
    
    
    T2H_T2
    saveas(fig1,'figures/plots/'+string(ii)+'_fig1_H_S_.png')
    m = [cTemp', T2];
    writematrix(m,'figures/data/'+string(ii)+'_H_S_.csv')
    close all hidden
%     clear
%     clc
    
end