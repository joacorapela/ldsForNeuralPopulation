function [pll_slc_fromAvctrl,pll_slc_fromModelPred,pll_slc_fromDataAv,...
    pll_cntrl_fromAvctrl,pll_cntrl_fromModelPred,pll_cntrl_fromDataAv,pll_cntrl_fromLowdim] = ...
    makePLLs(predictedNeurons,tSilencingLag,timewindow,TrY,TrYcntrl,TrYcntrl_red,Ypred_slc,Ypred,condi,...
    pll_slc_fromAvctrl,pll_slc_fromModelPred,pll_slc_fromDataAv,...
    pll_cntrl_fromAvctrl,pll_cntrl_fromModelPred,pll_cntrl_fromDataAv,pll_cntrl_fromLowdim);
% in comparisons to averges, the trial itseld is excluded, and others averaged
for NeuronNum = predictedNeurons % 54,3, animal 2
    
    pll_slc_fromAvctrl(tSilencingLag,NeuronNum) = 0; % use a randon control trial to predict it, average pll for all control trials
    pll_slc_fromModelPred(tSilencingLag,NeuronNum) = 0; % same trial prediction
    pll_slc_fromDataAv(tSilencingLag,NeuronNum) = 0; % from average across trials from data (silecing trials)
    
    pll_cntrl_fromAvctrl(tSilencingLag,NeuronNum) = 0;
    pll_cntrl_fromModelPred(tSilencingLag,NeuronNum) = 0; % same trial prediction
    pll_cntrl_fromDataAv(tSilencingLag,NeuronNum) = 0; % from average across trials from data (control trials) -- if 0 returns inf and nan, fix [TODO]
    pll_cntrl_fromLowdim(tSilencingLag,NeuronNum) = 0;
    % also, the trial itseld should be excluded, average others  (for both model control and av data) [TODO]
    
    
    N_slc_fromAvctrl = 0;
    N_slc_fromModelPred = 0;
    N_slc_fromDataAv = 0;
    
    N_cntrl_fromAvctrl = 0;
    N_cntrl_fromModelPred = 0;
    N_cntrl_fromDataAv = 0;
    N_cntrl_fromLowdim = 0;
    
    for slctr = 1:size(TrY(:,timewindow,NeuronNum),1)
        elsetr = setdiff(1:size(TrY(:,timewindow,NeuronNum),1),slctr);
        for timepoint = 1:numel(timewindow)
            
            x = TrY(slctr,timewindow(timepoint),NeuronNum);
            lambda_fromModelPred = Ypred_slc(slctr,timewindow(timepoint),NeuronNum);
            lambda_fromAvctrl = Ypred(:,timewindow(timepoint),NeuronNum); % this is a vec, not a number, for all trials
            lambda_fromDataAv = mean(TrY(elsetr,timewindow(timepoint),NeuronNum),1);
            
            if ~isnan(lambda_fromModelPred) % && lambda_fromModelPred < 1000% skips if lambda = 0, or if there is no prediction (nan) -- are they the same?
                if lambda_fromModelPred <= 0
                    keyboard
                    error('')
                end
                pll_slc_fromModelPred(tSilencingLag,NeuronNum) = nansum([pll_slc_fromModelPred(tSilencingLag,NeuronNum) , ...
                    x*log(lambda_fromModelPred)-lambda_fromModelPred-log(factorial(x))]);
                N_slc_fromModelPred = N_slc_fromModelPred + 1;
            end
            
            if any(~isnan(lambda_fromAvctrl)) % if any trial is fine (or shoud be all?)
                if any(lambda_fromAvctrl<= 0)
                    keyboard
                    error('')
                end
                pll_slc_fromAvctrl(tSilencingLag,NeuronNum) = nansum([pll_slc_fromAvctrl(tSilencingLag,NeuronNum) , ...
                    nanmean(x.*log(lambda_fromAvctrl)-lambda_fromAvctrl-log(factorial(x)))]);
                N_slc_fromAvctrl = N_slc_fromAvctrl +1;
            end
            
            if lambda_fromDataAv > 0 % this can't be nan, do we increase n inside or outside? probably outside?
                pll_slc_fromDataAv(tSilencingLag,NeuronNum) = nansum([pll_slc_fromDataAv(tSilencingLag,NeuronNum) , ...
                    x.*log(lambda_fromDataAv)-lambda_fromDataAv-log(factorial(x))]);
                N_slc_fromDataAv = N_slc_fromDataAv +1; % check
            end
        end
    end
    
    % average per number of trials and time points, to make
    % comparable
    
    pll_slc_fromModelPred(tSilencingLag,NeuronNum) = pll_slc_fromModelPred(tSilencingLag,NeuronNum)/N_slc_fromModelPred;%(numel(timewindow)*size(TrY(:,timewindow,NeuronNum),1));
    pll_slc_fromAvctrl(tSilencingLag,NeuronNum) = pll_slc_fromAvctrl(tSilencingLag,NeuronNum)/N_slc_fromAvctrl;%(numel(timewindow)*size(TrY(:,timewindow,NeuronNum),1));
    pll_slc_fromDataAv(tSilencingLag,NeuronNum) = pll_slc_fromDataAv(tSilencingLag,NeuronNum)/N_slc_fromDataAv;%(numel(timewindow)*size(TrY(:,timewindow,NeuronNum),1));
    
    for ctrtr = 1:size(TrYcntrl(:,timewindow,NeuronNum),1)
        elsetr = setdiff(1:size(TrYcntrl(:,timewindow,NeuronNum),1),ctrtr);
        for timepoint = 1:numel(timewindow)
            
            x = TrYcntrl(ctrtr,timewindow(timepoint),NeuronNum);
            lambda_fromModelPred = Ypred(ctrtr,timewindow(timepoint),NeuronNum);
            lambda_fromAvctrl = Ypred(elsetr,timewindow(timepoint),NeuronNum); % this is a vector
            lambda_fromDataAv = mean(TrYcntrl(elsetr,timewindow(timepoint),NeuronNum),1);
            lambda_fromLowdim = nanmean(TrYcntrl_red(elsetr,timewindow(timepoint),NeuronNum),1);
            
            clear cond_phrase
            if isnan(condi)
                cond_phrase = 1;
            elseif condi == 1
                cond_phrase = lambda_fromDataAv > 0 &&  ~isnan(lambda_fromModelPred) && lambda_fromLowdim > 0 && any(~isnan(lambda_fromAvctrl)); % same condition for multiple plots
            end
            
            if cond_phrase
                
                if ~isnan(lambda_fromModelPred)
                    if lambda_fromModelPred <= 0
                        keyboard
                        error('')
                    end
                    pll_cntrl_fromModelPred(tSilencingLag,NeuronNum) = nansum([pll_cntrl_fromModelPred(tSilencingLag,NeuronNum) , ...
                        x * log(lambda_fromModelPred)-lambda_fromModelPred-log(factorial(x))]);
                    N_cntrl_fromModelPred = N_cntrl_fromModelPred + 1;
                end
                
                
                if any(~isnan(lambda_fromAvctrl)) % if any trial is fine (or shoud be all?)
                    if any(lambda_fromAvctrl<= 0)
                        keyboard
                        error('')
                    end
                    pll_cntrl_fromAvctrl(tSilencingLag,NeuronNum) = nansum([pll_cntrl_fromAvctrl(tSilencingLag,NeuronNum) , ...
                        nanmean( x.*log(lambda_fromAvctrl)-lambda_fromAvctrl-log(factorial(x)) )]);
                    N_cntrl_fromAvctrl = N_cntrl_fromAvctrl + 1;
                end
                
                
                if lambda_fromDataAv > 0
                    pll_cntrl_fromDataAv(tSilencingLag,NeuronNum) = nansum([pll_cntrl_fromDataAv(tSilencingLag,NeuronNum) , ...
                        x.*log(lambda_fromDataAv)-lambda_fromDataAv-log(factorial(x))]);
                    N_cntrl_fromDataAv = N_cntrl_fromDataAv + 1;
                end
                
                
                if  lambda_fromLowdim > 0
                    pll_cntrl_fromLowdim(tSilencingLag,NeuronNum) = nansum([pll_cntrl_fromLowdim(tSilencingLag,NeuronNum) , ...
                        x.*log(lambda_fromLowdim)-lambda_fromLowdim-log(factorial(x))]);
                    N_cntrl_fromLowdim = N_cntrl_fromLowdim + 1;
                end
                
                
            end
        end
    end
    pll_cntrl_fromModelPred(tSilencingLag,NeuronNum) = pll_cntrl_fromModelPred(tSilencingLag,NeuronNum)/N_cntrl_fromModelPred;%(numel(timewindow)*size(TrYcntrl(:,timewindow,NeuronNum),1));
    pll_cntrl_fromDataAv(tSilencingLag,NeuronNum) =  pll_cntrl_fromDataAv(tSilencingLag,NeuronNum)/N_cntrl_fromDataAv;%(numel(timewindow)*size(TrYcntrl(:,timewindow,NeuronNum),1));
    pll_cntrl_fromAvctrl(tSilencingLag,NeuronNum) = pll_cntrl_fromAvctrl(tSilencingLag,NeuronNum)/N_cntrl_fromAvctrl;%(numel(timewindow)*size(TrYcntrl(:,timewindow,NeuronNum),1));
    pll_cntrl_fromLowdim(tSilencingLag,NeuronNum) = pll_cntrl_fromLowdim(tSilencingLag,NeuronNum)/N_cntrl_fromLowdim;%(numel(timewindow)*size(TrYcntrl(:,timewindow,NeuronNum),1));
end