% Decoding the stimulus identity from V1 activity in the presence and
% absence of feedback, using a nonlinear (quadratic) svm.

% an80, (much)better when cleaned. 60ms looks kinda the same
%  or pl5- take 80-120
%
% *** IMPORTANT INFO
% Results depend on whether or not the cells are cleaned and removed.
% Altough the patterns are similar (even significance), and seems that both
% removing and cleaning looks the best, after that only remove, after that
% no remove or clean. (can take both remove and clean to be consistent with
% cormatrix_prediction and timeconstant analysis.)
% Also note that cleaning would be done with 80ms bins, which is not the
% same as other analysis with 150ms. Default is both clean and remove for now.
% it also depends on nfolds. Default 2
% Also, should there been only 1 divide or many? For multiple re-draws, rep should be set
% to >1. But this makes the error bars very small. (Does it make sense?),
% but eliminates sensitivity to random generator seed
%

%% CV with folds, and add time axis (steps of mean and ste). Add pvalues
% load 70ms, clean and remove, 
% later fix the ci and the tests: can do for individual neurons? or Leva one trial out and
% individual trials?

% For his code, the size of analysis bin matters
clc
nfold=3; % 3 works best in showing the difference -  2 is the best
numanimals=6;
numtimepoints = 8;
nrep = 20;
Mbs= nan(1,numtimepoints);
Mls= nan(1,numtimepoints);
Ebs= nan(1,numtimepoints);
Els= nan(1,numtimepoints);

for timepoint = 1:numtimepoints
    
    
    ClAc_bs = nan(numanimals,nfold);
    ClAc_ls = nan(numanimals,nfold);
    for animalnum=1:numanimals
        for rep=1:nrep
            
            targetcell = V1cells(find(cellfun(@(x) x.simulcode,V1cells) == animalnum ));
            % ntrials*ncells:
            gobs = cell2mat(cellfun(@(x) x.laAbs.go{timepoint},targetcell,'UniformOutput',0));
            nogobs = cell2mat(cellfun(@(x) x.laAbs.nogo{timepoint},targetcell,'UniformOutput',0));
            gols = cell2mat(cellfun(@(x) x.laAls.go{timepoint},targetcell,'UniformOutput',0));
            nogols = cell2mat(cellfun(@(x) x.laAls.nogo{timepoint},targetcell,'UniformOutput',0));
            %     % for this, load plot 5ms
            %     gobs = cell2mat(cellfun(@(x) sum(x.nbs.go(:,200+(16:24)),2),targetcell,'UniformOutput',0));
            %     nogobs = cell2mat(cellfun(@(x) sum(x.nbs.nogo(:,200+(16:24)),2),targetcell,'UniformOutput',0));
            %     gols = cell2mat(cellfun(@(x) sum(x.nls.go{timpepoint}(:,200+(16:24)),2),targetcell,'UniformOutput',0));
            %     nogols = cell2mat(cellfun(@(x) sum(x.nls.nogo{timpepoint}(:,200+(16:24)),2),targetcell,'UniformOutput',0));
            
            bsall = [gobs;nogobs];
            lsall = [gols;nogols];
            
            
            % equalize trials
            gobs(:,isnan(sum(bsall,1))) = [];
            nogobs(:,isnan(sum(bsall,1))) = [];
            gols(:,isnan(sum(lsall,1))) = [];
            nogols(:,isnan(sum(lsall,1))) = [];
            
            mintrialnum = min([size(gobs,1),size(nogobs,1),size(gols,1),size(nogols,1)]);
            gobs = gobs(1:mintrialnum,:);
            nogobs = nogobs(1:mintrialnum,:);
            gols = gols(1:mintrialnum,:);
            nogols = nogols(1:mintrialnum,:);
            
            
            rng(rep)
           % rng(3)
            % LOO might be better
            cvind=crossvalind('Kfold',size(gobs,1),nfold);
            %cvind=crossvalind('LeaveMOut',size(gobs,1),1);
            
            for fold=1:nfold
                traininds = (find(cvind ~= fold));
                testinds = (find(cvind == fold));
                % baseline
                mdl_bs = fitcsvm([gobs(traininds,:);nogobs(traininds,:)],...
                    [zeros(numel(traininds),1);ones(numel(traininds),1)],'KernelFunction',...
                    'polynomial','PolynomialOrder',2,'Prior','uniform','Standardize',0);
%                 mdl_bs = fitcdiscr([gobs(traininds,:);nogobs(traininds,:)],...
%                     [zeros(numel(traininds),1);ones(numel(traininds),1)],'DiscrimType',...
%                     'diagQuadratic','Prior','empirical');
                predicted_label = predict(mdl_bs,[gobs(testinds,:);nogobs(testinds,:)]);
                real_label = [zeros(numel(testinds),1);ones(numel(testinds),1)];
                % check classification accuracy
                % ClAc_bs(end+1) = sum(~xor(predicted_label,real_label))/numel(predicted_label);
                ClAc_bs(animalnum,fold,rep) = sum(~xor(predicted_label,real_label))/numel(predicted_label);
                
                % laser
                mdl_ls = fitcsvm([gols(traininds,:);nogols(traininds,:)],...
                    [zeros(numel(traininds),1);ones(numel(traininds),1)],'KernelFunction',...
                    'polynomial','PolynomialOrder',2,'Prior','uniform','Standardize',0);
%                 mdl_ls = fitcdiscr([gols(traininds,:);nogols(traininds,:)],...
%                     [zeros(numel(traininds),1);ones(numel(traininds),1)],'DiscrimType',...
%                     'diagQuadratic','Prior','empirical');
                predicted_label = predict(mdl_ls,[gols(testinds,:);nogols(testinds,:)]);
                real_label = [zeros(numel(testinds),1);ones(numel(testinds),1)];
                % check classification accuracy
%                 ClAc_ls(end+1) = sum(~xor(predicted_label,real_label))/numel(predicted_label);
                ClAc_ls(animalnum,fold,rep) = sum(~xor(predicted_label,real_label))/numel(predicted_label);
            end
            
        end
    end
   ClAc_bs = squeeze(nanmean(ClAc_bs,2));
   ClAc_ls = squeeze(nanmean(ClAc_ls,2));
    
%        ClAc_bs = nanmean(ClAc_bs,2);
%        ClAc_ls =nanmean(ClAc_ls,2);
    
    signrank(reshape(ClAc_bs,1,[]),reshape(ClAc_ls,1,[]))
    
    Mbs(timepoint) = nanmean(reshape(ClAc_bs,1,[]));
    Mls(timepoint) = nanmean(reshape(ClAc_ls,1,[]));
    
    Ebs(timepoint) = nanstd(reshape(ClAc_bs,1,[]))/sqrt(numel(ClAc_bs));
    Els(timepoint) = nanstd(reshape(ClAc_ls,1,[]))/sqrt(numel(ClAc_ls));
    
end
figure;shadedErrorBar([],Mbs,Ebs,'k',1);
hold on;
shadedErrorBar([],Mls,Els,'b',1)

wdt = 0.2;
figure;
pb=bar((1:numtimepoints)-wdt/2,Mbs,wdt,'FaceColor','k','EdgeColor','none'); hold on
errorbar((1:numtimepoints)-wdt/2,Mbs,Ebs,'k.'); hold on
pl=bar((1:numtimepoints)+wdt/2,Mls,wdt,'FaceColor','b','EdgeColor','none'); hold on
errorbar((1:numtimepoints)+wdt/2,Mls,Els,'b.')

ylim([0.5,1])
ax1 = gca;
ax1.XTickLabel=[];
ylabel('stimulus decoding accurancy from V1 cells (nonlinear SVM)');

xlim([0.5 3.5])
xlim([0 10])
ax1.XTick = [1 2 3 4 5 6 7 8];
%ax1.XTickLabel={'0ms-80ms','60ms-140ms','120ms-200ms'};
ax1.XTickLabel={'T1','T2','T3','T4','T5','T6','T7','T8'};
legend([pb,pl],{'baseline','laser'});
