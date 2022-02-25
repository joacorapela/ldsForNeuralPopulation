% % check predctions _ no cross validation 
% 
% 
% d = load(fullfile(config.summarymatfile.folder,config.summarymatfile.name))
% 

figure
% need to separate silencing from non-silencing trials
for i = 1:3
    if i ==1
        %fname = ('JointSep_lag0_PLDSfitRes_22_02_15_20_15_17_RND0_onlyCorrect_exGrooming.mat')
        fname = 'V1_PLDSfitRes_22_02_15_14_12_56_DiscInput_onlyCorrect_exGrooming';
    elseif i == 2
       
        fname = 'V1_PLDSfitRes_22_02_15_14_46_08_DiscInput_onlyCorrect_exGrooming_RND2.mat';
        fname = 'V1_PLDSfitRes_22_02_15_16_32_26_DiscInput_onlyCorrect_exGrooming_RND3.mat';
        fname = 'V1_PLDSfitRes_22_02_15_18_13_37_DiscInput_onlyCorrect_exGrooming_RND4.mat';
        
    elseif i == 3
        fname = 'V1_PLDSfitRes_22_02_15_14_54_17_DiscInput_onlyCorrect_exGrooming_RND9_14states';
    end
load(fname)

heldoutN = 30
pred_bs = nan(length(seq),seq(1).T); % trials*timebins
yOrig_bs = nan(length(seq),seq(1).T); % trials*timebins
pred_ls = nan(length(seq),seq(1).T); % trials*timebins
yOrig_ls = nan(length(seq),seq(1).T); % trials*timebins

for tr = 1:length(seq)
   
    % for trial tr:
    z = params.model.C(heldoutN,:) *seq(tr).posterior.xsm + ...
        params.model.d(heldoutN);
    %if sum(sum(seq(tr).u(3:end,:),1))>0
        if sum(sum(seq(tr).u(2,:),1)) >0
    pred_ls(tr,:) = exp(z);
    yOrig_ls(tr,:) = seq(tr).y(heldoutN,:);
       % end
    else
    pred_bs(tr,:) = exp(z);
    yOrig_bs(tr,:) = seq(tr).y(heldoutN,:);
    end
end
% pred is actually the expectation of the poisson distribution
% gamma is hard coded here

s = subplot(3,2,1+2*(i-1));plot(nanmean(pred_bs,1));hold on;plot(nanmean(yOrig_bs,1));legend({'bs prediction','bs original trace'})
s.Title.String = fname;s.Title.Interpreter = 'none';s.TitleFontWeight='normal';
 subplot(3,2,2+2*(i-1));plot(nanmean(pred_ls,1));hold on;plot(nanmean(yOrig_ls,1));legend({'ls prediction','ls original trace'})
end