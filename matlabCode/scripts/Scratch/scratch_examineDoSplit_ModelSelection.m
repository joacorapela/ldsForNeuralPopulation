
% local results, with laplace inference
nosplit = load(fullfile('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/VL61/trial_based_LONO/50msBins',...
    'ModelSelection-Fold1_V1_PLDSfitRes_21_09_02_19_15_59.mat'));
split = load(fullfile('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/VL61/trial_based_LONO/50msBins',...
    'ModelSelection-Fold1_V1_PLDSfitRes_21_09_02_19_18_53.mat'));

%%
ll = nan(1,50);
for st = 1:numel(ll)
    try
        if isfield(split.Allmodels{st},'Avtrial_ll')
            ll(st) = split.Allmodels{st}.Avtrial_ll;
        end
    catch
    end
end

figure;plot(1:50,ll,'b.-');
ll = nan(1,50);
for st = 1:numel(ll)
    try
        if isfield(nosplit.Allmodels{st},'Avtrial_ll')
            ll(st) = nosplit.Allmodels{st}.Avtrial_ll;
        end
    catch
    end
end

hold on;plot(1:50,ll,'r.-');
%%
vb = nan(1,50);
for st = 1:numel(vb)
    try
        vbs = split.Allmodels{st}.varBound;
        vbs = vbs(~isnan(vbs));
        vb(st) = vbs(end);
    catch
    end
end

figure;plot(1:50,vb,'b.-');

vb = nan(1,50);
for st = 1:numel(vb)
    try
        vbs = nosplit.Allmodels{st}.varBound;
        vbs = vbs(~isnan(vbs));
        vb(st) = vbs(end);
    catch
    end
end

hold on;plot(1:50,vb,'r.-');
legend('split','nosplit')