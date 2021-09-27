% Axis lims for MPV17, 25 ms

res = load('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/MPV17/trial_based_MS_split_1/25msBins/Fold1_V1_PLDSfitRes_21_09_16_18_09_05.mat');
tr_list_bs = [];
for tr = 1:length(res.seq)
  %  if numel(find(res.seq(tr).u(4,:)==0))
    if ~numel(find(res.seq(tr).u(3:end,:)==0)) & numel(find(res.seq(tr).u(1,:)==1))
        tr_list_bs = [tr_list_bs,tr];
    end
end

tr_list_ls = [];
for tr = 1:length(res.seq)
    if numel(find(res.seq(tr).u(4,:)==0))
        tr_list_ls = [tr_list_ls,tr];
    end
end

fb_av_bs = nanmean(cell2mat(arrayfun(@(x) x.u(4,:),res.seq(tr_list_bs),'UniformOutput',0)'),1);
stim_av_bs = nanmean(cell2mat(arrayfun(@(x) x.u(1,:),res.seq(tr_list_bs),'UniformOutput',0)'),1);
fb_av_ls = nanmean(cell2mat(arrayfun(@(x) x.u(4,:),res.seq(tr_list_ls),'UniformOutput',0)'),1);
stim_av_ls = nanmean(cell2mat(arrayfun(@(x) x.u(1,:),res.seq(tr_list_ls),'UniformOutput',0)'),1);

y_av_bs = nan(size(res.seq(1).y));
for neuron = 1:size(res.seq(1).y,1)
y_av_bs(neuron,:) = nanmean(cell2mat(arrayfun(@(x) x.y(neuron,:),res.seq(tr_list_bs),'UniformOutput',0)'),1);
end

x_av_bs = nan(size(res.seq(1).posterior.xsm));
x_ci_bs = nan(size(res.seq(1).posterior.xsm));
for st = 1:size(res.seq(1).posterior.xsm,1)
x_av_bs(st,:) = nanmean(cell2mat(arrayfun(@(x) x.posterior.xsm(st,:),res.seq(tr_list_bs),'UniformOutput',0)'),1);
temp = cell2mat(arrayfun(@(x) x.posterior.xsm(st,:),res.seq(tr_list_bs),'UniformOutput',0)');
x_ci_bs(st,:) = 2*nanstd(temp)/sqrt(size(temp,1));
end

y_av_ls = nan(size(res.seq(1).y));
for neuron = 1:size(res.seq(1).y,1)
y_av_ls(neuron,:) = nanmean(cell2mat(arrayfun(@(x) x.y(neuron,:),res.seq(tr_list_ls),'UniformOutput',0)'),1);
end

x_av_ls = nan(size(res.seq(1).posterior.xsm));
x_ci_ls = nan(size(res.seq(1).posterior.xsm));
for st = 1:size(res.seq(1).posterior.xsm,1)
x_av_ls(st,:) = nanmean(cell2mat(arrayfun(@(x) x.posterior.xsm(st,:),res.seq(tr_list_ls),'UniformOutput',0)'),1);
temp = cell2mat(arrayfun(@(x) x.posterior.xsm(st,:),res.seq(tr_list_ls),'UniformOutput',0)');
x_ci_ls(st,:) = 2*nanstd(temp)/sqrt(size(temp,1));
end


figure;plot(fb_av_bs);
hold on; plot(stim_av_bs);
hold on; plot(y_av_bs','k')

figure;
s_ind = find(stim_av_bs);
l_ind = find(~fb_av_ls);
for st = 1:8
    s = subplot(1,8,st);
%     plot(fb_av_bs);
%     hold on; plot(fb_av_ls);
%     hold on; plot(stim_av_bs);
%     hold on; plot(stim_av_ls);
    hold on;patch([s_ind(1),s_ind(end),s_ind(end),s_ind(1)],[-1 -1,1,1],'y','FaceAlpha',0.1,...
    'EdgeAlpha',0)
    hold on;patch([l_ind(1),l_ind(end),l_ind(end),l_ind(1)],[-1 -1,1,1],'r','FaceAlpha',0.1,...
    'EdgeAlpha',0)
    hold on; plot(x_av_bs(st,:),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
    hold on; patch([1:80,80:-1:1],[x_av_bs(st,:)-x_ci_bs(st,:),fliplr(x_av_bs(st,:)+x_ci_bs(st,:))],[0.5,0.5,0.5],...
        'FaceAlpha',0.2,'EdgeAlpha',0)
    hold on; plot(x_av_ls(st,:),'Color',[0,0.5,1],'LineWidth',1.5)
    hold on; patch([1:80,80:-1:1],[x_av_ls(st,:)-x_ci_ls(st,:),fliplr(x_av_ls(st,:)+x_ci_ls(st,:))],[0,0.5,1],...
        'FaceAlpha',0.2,'EdgeAlpha',0)
    ylim([-0.03 0.04])    
    s.Title.String = ['state ',num2str(st)];
    s.Title.FontWeight = 'normal';
    
    %%% This is the onyl part that is bin related. comment out if changing
    %%% bin size, default 25 ms
    xlim([30 70])
    s.XTick = [40 60];
    s.XTickLabel = [0 ,500];
    xlabel('time (ms)')
    
    if st>1
        s.YTickLabel = [];
    end
end
