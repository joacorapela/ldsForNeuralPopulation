
cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/VL63/trial_based/50msBins')

area = 'LM';
nlags = 11; % delays (including zero)
ninputs = 4;
do_plot = 1;
nexamples = 10; % number of example cells to plot

PLDSresFiles = dir([area,'*.mat']);
% by default takes the latest file
if length(PLDSresFiles)>1
    PLDSresFiles = PLDSresFiles(1);
end
res = load(PLDSresFiles.name);
RF = cell(1,size(res.params.model.C,1)); % number of neurons
for n = 1:length(RF)
    RF{n} = nan(nlags,ninputs);
end
RF = calcERFforWindow(res,RF,nlags,ninputs);


if do_plot
    for n = 1:nexamples
        figure;
        hold on; plot(RF{n}(:,1),'Color',[0 0.7 0 ]);% vis go
        hold on; plot(RF{n}(:,2),'Color',[0.7 0 0 ]);% vis nogo
        hold on; plot(RF{n}(:,3),'g');% laser go
        hold on; plot(RF{n}(:,4),'r');% laser nogo
        legend({'vg','vn','lg','ln'})
        
    end
end

function RF = calcERFforWindow(res,RF,nlags,ninputs)
for n = 1:size(RF,2) % n is number of neurons
    for in = 1:ninputs
        for i=1:nlags
            RF{n}(i,in) = res.params.model.C(n,:)*((res.params.model.A)^(i-1))*res.params.model.B(:,in);
        end
    end
end
end

% figure;
% h = plot(params.model.B);
% set(h, {'color'}, {[0 0.7 0 ]; [0.7 0 0 ]; [0 1 0]; [1 0 0]});
% legend({'vg','vn','lg','ln'})


