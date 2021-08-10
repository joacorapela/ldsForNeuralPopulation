

% fit on training dataset
LONOpath = '/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/VL61/trial_based_LONO/50msBins';
LONOname_V1 = 'Fold1_V1_PLDSfitRes_21_08_09_14_27_24.mat';
LONOname_LM = 'Fold1_LM_PLDSfitRes_21_08_09_15_19_33.mat';

% fit on fulldataset
Fullpath = '/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/VL61/trial_based/50msBins';
Fullname_V1 = 'V1_PLDSfitRes_21_07_23_18_52_55.mat';
Fullname_LM = 'LM_PLDSfitRes_21_07_23_20_00_55.mat';

resLONO = load(fullfile(LONOpath,LONOname_V1));
resFull = load(fullfile(Fullpath,Fullname_V1));

figure;subplot(1,2,1);imagesc(resLONO.params.model.A);
subplot(1,2,2);imagesc(resFull.params.model.A);

figure;subplot(1,2,1);imagesc(resLONO.params.model.B);
subplot(1,2,2);imagesc(resFull.params.model.B);

figure;plot(cell2mat(arrayfun(@(x) x.u(1,:) - x.u(2,:), resLONO.seq,'UniformOutput',0)')')

% checking trial averaged state variables:
for state=1:9
    figure
    hold on;plot(mean(cell2mat(arrayfun(@(x) x.posterior.xsm(state,:), resLONO.seq,'UniformOutput',0)'),1),'r')
    hold on;plot(mean(cell2mat(arrayfun(@(x) x.posterior.xsm(state,:), resFull.seq,'UniformOutput',0)'),1),'k')
end

figure;subplot(1,2,1);imagesc(resLONO.params.model.C);
subplot(1,2,2);imagesc(resFull.params.model.C);
%% receptive field of neurons calculated from the full or LONO model: very similar

n=11;   %neuron number 
RF = [];

for in = 1:4
    for i=1:10
        RF{n}(i,in) = resFull.params.model.C(n,:)*((resFull.params.model.A)^(i-1))*resFull.params.model.B(:,in);
    end
end


figure;plot(RF{n}(:,1),'g');hold on;plot(RF{n}(:,2),'r');
hold on;plot(RF{n}(:,3),'g');hold on;plot(RF{n}(:,4),'r');

RF = [];

for in = 1:4
    for i=1:10
        RF{n}(i,in) = resLONO.params.model.C(n,:)*((resLONO.params.model.A)^(i-1))*resLONO.params.model.B(:,in);
    end
end



hold on;plot(RF{n}(:,1),'g--');hold on;plot(RF{n}(:,2),'r--');
hold on;plot(RF{n}(:,3),'g--');hold on;plot(RF{n}(:,4),'r--');