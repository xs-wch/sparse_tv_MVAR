%%%% This version to run the test program on cluster
clear all
close all
ntest = 100;
A_all = cell(ntest ,1);
stv_mvar_all = cell(ntest ,1);
SNR_test = [20,15,10];
changeamp = [0.7,0.5,0.3,0.1];
matlabpool open
for iSNR = 1:length(SNR_test)
    for iamp = 1:length(changeamp)
        parfor itest = 1:ntest
            stv_mvar = sparse_time_varying_MVAR_EMpara(3,SNR_test(iSNR),changeamp(iamp),[]);
            A = stv_mvar.estimate_model_SAEM();
            
%             figure
%             for i = 1:6
%                 for j = 1:6
%                     subplot(6,6,(i-1)*6+j)
%                     chanindex = [i,j];
%                     plot(1:499, squeeze(A(chanindex(1),chanindex(2),:)))
%                     hold on
%                     plot([1:500], [repmat(stv_mvar.synthetic_A1(chanindex(1),chanindex(2)),1,200), repmat(stv_mvar.synthetic_A2(chanindex(1),chanindex(2)),1,300)],'r-')
%                     ylim([-1, 1])
%                     xlim([0, 500])
%                     xlabel('time(ms)')
%                     
%                 end
%             end
            
            A_all{itest} = A;
            stv_mvar_all{itest} = stv_mvar;

        end
        save(['stv_mvar_SNR',num2str(SNR_test(iSNR)),'_amp',num2str(changeamp(iamp)),'.mat'],'A_all','stv_mvar_all');
    end
end
matlabpool close