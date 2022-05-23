
animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};

% orthogonalization not relevant here
% This code can chooose between leave-one-area-out, simialr to analyse_predict_Xtrials
% andleave-one-neuron, as in analyse_predict_Xtrials_LONO, and can
% when slc trials used, both input is received and y of LM is observed
% (which is silenced)
savedir = ('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/LONO_LOAO');
cd(savedir)
diary(datestr(now))
diary on

rng(1,'twister');
CntrlOnly = 1 % 0 or not
nfold = 5;
leaveout = 'neuron' % neuron or area: all V1


if CntrlOnly
    UInd = 1;
else
    UInd = 1:9;
end

for animali = 1:6 
    
    disp(['animal = ',num2str(animali)]);
        
    if CntrlOnly
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext_CntrlOnly/17msBins/'])
        filelist = dir(['*splt1_RND0_onlyCorrect_exGrooming_go_',num2str(nfold),'FoldXval*.mat']);
    else
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
        filelist = dir(['*splt1_onlyCorrect_exGrooming_go_',num2str(nfold),'FoldXval*.mat']);
    end
    res = load(filelist(end).name);
    
    NV1Cells = find(sum(res.Xval.params{1}.model.C(:,1:8),2)==0,1)-1%27;% all 37 % make dynamic maybe try animals with more lm cells
    pvcell =[]%36
    
    Lag = cell(1,8);
    Ypred = cell(1,nfold);
    TrYcntrl = cell(1,nfold);
    for i =1:8
        Lag{i}.Ypred_slc = cell(1,nfold);
        Lag{i}.TrY = cell(1,nfold);
    end
    % maybe push tSilencinglag inside to reduce 7 extra control
    % calculation, an save the 4 traces afterwards. Then load for the
    % next part.
    for fold = 1:nfold
        if strcmp(leaveout, 'area')
            warning off
            [Ypred{fold},TrYcntrl{fold}] = makeYs_cntrl(res,NV1Cells,pvcell,UInd,fold);
            for tSilencingLag = 1:8
                [Lag{tSilencingLag}.Ypred_slc{fold},Lag{tSilencingLag}.TrY{fold}] = makeYs_slc(res,NV1Cells,pvcell,UInd,fold,tSilencingLag);
            end
            % [Ypred{fold},Ypred_slc{fold},TrYcntrl{fold},TrY{fold}] = makeYs(res,NV1Cells,pvcell,UInd,fold,tSilencingLag);
        elseif strcmp(leaveout, 'neuron')
            warning off % warning of sigular matrix in laplaceinferencecore;
            [Ypred{fold},TrYcntrl{fold}] = makeYs_LONO_cntrl(res,UInd,fold);
            for tSilencingLag = 1:8
                [Lag{tSilencingLag}.Ypred_slc{fold},Lag{tSilencingLag}.TrY{fold}] = makeYs_LONO_slc(res,UInd,fold,tSilencingLag);
            end
            %[Ypred{fold},Ypred_slc{fold},TrYcntrl{fold},TrY{fold}] = makeYs_LONO(res,UInd,fold,tSilencingLag);
        end
    end
    %%% here to the number of folds
    
    Ypred = cat(1,Ypred{1},Ypred{2},Ypred{3},Ypred{4},Ypred{5});
    TrYcntrl = cat(1,TrYcntrl{1},TrYcntrl{2},TrYcntrl{3},TrYcntrl{4},TrYcntrl{5});
    for tSilencingLag = 1:8
        Lag{tSilencingLag}.Ypred_slc = cat(1, Lag{tSilencingLag}.Ypred_slc{1}, Lag{tSilencingLag}.Ypred_slc{2}, Lag{tSilencingLag}.Ypred_slc{3}, Lag{tSilencingLag}.Ypred_slc{4}, Lag{tSilencingLag}.Ypred_slc{5});
        Lag{tSilencingLag}.TrY = cat(1, Lag{tSilencingLag}.TrY{1}, Lag{tSilencingLag}.TrY{2}, Lag{tSilencingLag}.TrY{3}, Lag{tSilencingLag}.TrY{4}, Lag{tSilencingLag}.TrY{5});
    end
    
    %%%% saving
    savename  = ['animal',num2str(animali),'_',leaveout,'_','CntrlOnly',num2str(CntrlOnly)];
    save(fullfile(savedir,savename),'Ypred','TrYcntrl','Lag','NV1Cells','pvcell','nfold')
end

diary off
