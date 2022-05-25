function [dhat,dtrue,dhat_chance,dtrue_chance,dhat_boot,dtrue_boot] = ...
    makeDsX(dhat,dtrue,dhat_chance,dtrue_chance,dhat_boot,dtrue_boot,nboot,timewindow,NV1Cells,Ypred,Ypred_slc,TrYcntrl,TrY,dMethod,activity_thresh,tSilencingLag)
if strcmp(dMethod,'lda')
    error('not implemented')
    %     X = [squeeze(nanmean(Ypred(:,timewindow,:),2));squeeze(nanmean(Ypred_slc(:,timewindow,:),2))];
    %     Y = [zeros(size(Ypred,1),1);ones(size(Ypred_slc,1),1)];
    %     mdl = fitcdiscr(X,Y,'discrimType','pseudolinear');
    %     dhat(tSilencingLag,:) = mdl.Coeffs(1,2).Linear;
    %
    %     X = [squeeze(nanmean(TrYcntrl(:,timewindow,:),2));squeeze(nanmean(TrY(:,timewindow,:),2))];
    %     Y = [zeros(size(TrYcntrl,1),1);ones(size(TrY,1),1)];
    %     mdl = fitcdiscr(X,Y,'discrimType','pseudolinear');
    %     dtrue(tSilencingLag,:) = mdl.Coeffs(1,2).Linear;
    %
    %
    %     dhat_chance = [];
    %     dtrue_chance = [];
    %
    %     dhat_boot = [];
    %     dtrue_boot = [];
    
elseif strcmp(dMethod,'mean')
    for NeuronNum = 1:NV1Cells
        
        if isnan(activity_thresh)
            %             dhat(tSilencingLag,NeuronNum) = nanmean(nanmean(Ypred_slc(:,timewindow,(NeuronNum)),2)) - nanmean(nanmean(Ypred(:,timewindow,(NeuronNum)),2));
            %             dtrue(tSilencingLag,NeuronNum) = nanmean(nanmean(TrY(:,timewindow,(NeuronNum)),2)) - nanmean(nanmean(TrYcntrl(:,timewindow,(NeuronNum)),2));
            
            nCntrlTr = size(Ypred(:,timewindow,NeuronNum),1);
            for b = 1:nboot
                reord = randperm(nCntrlTr);
                half1 = reord(1:floor(numel(reord)/2));
                half2 = setdiff(reord,half1);
                
                dhat.p1(tSilencingLag,NeuronNum,b) = nanmean(nanmean(Ypred_slc(:,timewindow,(NeuronNum)),2)) - nanmean(nanmean(Ypred(half1,timewindow,(NeuronNum)),2));
                dhat.p2(tSilencingLag,NeuronNum,b) = nanmean(nanmean(Ypred_slc(:,timewindow,(NeuronNum)),2)) - nanmean(nanmean(Ypred(half2,timewindow,(NeuronNum)),2));
                dtrue.p1(tSilencingLag,NeuronNum,b) = nanmean(nanmean(TrY(:,timewindow,(NeuronNum)),2)) - nanmean(nanmean(TrYcntrl(half1,timewindow,(NeuronNum)),2));
                dtrue.p2(tSilencingLag,NeuronNum,b) = nanmean(nanmean(TrY(:,timewindow,(NeuronNum)),2)) - nanmean(nanmean(TrYcntrl(half2,timewindow,(NeuronNum)),2));
                
                dtrue_chance(tSilencingLag,NeuronNum,b) = nanmean(nanmean(TrYcntrl(half1,timewindow,(NeuronNum)),2)) - nanmean(nanmean(TrYcntrl(half2,timewindow,(NeuronNum)),2));
                dhat_chance(tSilencingLag,NeuronNum,b) = nanmean(nanmean(Ypred(half1,timewindow,(NeuronNum)),2)) - nanmean(nanmean(Ypred(half2,timewindow,(NeuronNum)),2));
            end
        else
            error('not implemented')
            %             if nanmean(nanmean(TrY(:,timewindow,(NeuronNum)),2)) > activity_thresh
            %                 dhat(tSilencingLag,NeuronNum) = nanmean(nanmean(Ypred_slc(:,timewindow,(NeuronNum)),2)) - nanmean(nanmean(Ypred(:,timewindow,(NeuronNum)),2));
            %                 dtrue(tSilencingLag,NeuronNum) = nanmean(nanmean(TrY(:,timewindow,(NeuronNum)),2)) - nanmean(nanmean(TrYcntrl(:,timewindow,(NeuronNum)),2));
            %             else
            %                 dhat(tSilencingLag,NeuronNum) = 0;
            %                 dtrue(tSilencingLag,NeuronNum) = 0;
            %             end
        end
        
        
        %         for b = 1:nboot
        %             slc_trials = datasample(1:numel(nanmean(Ypred_slc(:,timewindow,(NeuronNum)),2)),numel(nanmean(Ypred_slc(:,timewindow,(NeuronNum)),2)));
        %             cntrl_trials = datasample(1:numel(nanmean(Ypred(:,timewindow,(NeuronNum)),2)),numel(nanmean(Ypred(:,timewindow,(NeuronNum)),2)));
        %             cntrl_trials2 = datasample(1:numel(nanmean(Ypred(:,timewindow,(NeuronNum)),2)),numel(nanmean(Ypred(:,timewindow,(NeuronNum)),2)));
        %
        %             dhat_chance(tSilencingLag,NeuronNum,b) = nanmean(nanmean(Ypred(cntrl_trials,timewindow,(NeuronNum)),2)) - nanmean(nanmean(Ypred(cntrl_trials2,timewindow,(NeuronNum)),2)); % cntrl_trials2 or : for second term
        %             dtrue_chance(tSilencingLag,NeuronNum,b) = nanmean(nanmean(TrYcntrl(cntrl_trials,timewindow,(NeuronNum)),2)) - nanmean(nanmean(TrYcntrl(cntrl_trials2,timewindow,(NeuronNum)),2)); % cntrl_trials2 or : for second term
        %
        %             dhat_boot(tSilencingLag,NeuronNum,b) = nanmean(nanmean(Ypred_slc(slc_trials,timewindow,(NeuronNum)),2)) - nanmean(nanmean(Ypred(cntrl_trials,timewindow,(NeuronNum)),2));
        %             dtrue_boot(tSilencingLag,NeuronNum,b) = nanmean(nanmean(TrY(slc_trials,timewindow,(NeuronNum)),2)) - nanmean(nanmean(TrYcntrl(cntrl_trials,timewindow,(NeuronNum)),2));
        %         end
        
    end
end
end