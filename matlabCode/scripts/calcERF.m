
% run from results/animalName/binsize/tw folder

area = 'V1';
nlags = 11; % delays (including zero)
ninputs = 4;
do_plot = 1;
nexamples = 10; % number of example cells to plot

PLDSresFiles = dir([area,'*.mat']);
nreps = length(PLDSresFiles);

res = load(PLDSresFiles(1).name);
RF = cell(1,size(res.params.model.C,1)); % number of neurons
for n = 1:length(RF)
    RF{n} = nan(nlags,ninputs,nreps);
end

for rep = 1:nreps
    res = load(PLDSresFiles(rep).name);
    % uncomment to see the likelihood of each model
%    res.seq.posterior.varBound
    RF = calcERFforWindow(res,RF,nlags,ninputs,rep);
end

if do_plot
    for n = 1:nexamples
        figure;
        svg = subplot(3,2,1);
        vg = plot(squeeze(RF{n}(:,1,:)));% vis go
        set(vg, {'color'}, mat2cell(copper(nreps),ones(1,nreps)));
        ylim([min(min(squeeze(min(RF{n})))),max(max(squeeze(max(RF{n}))))]);
        svg.Title.String = 'vg';

        svng = subplot(3,2,2); 
        vng = plot(squeeze(RF{n}(:,2,:)));% vis nogo
        set(vng, {'color'}, mat2cell(copper(nreps),ones(1,nreps)));
        ylim([min(min(squeeze(min(RF{n})))),max(max(squeeze(max(RF{n}))))]);
        svng.Title.String = 'vng';
        
        slg = subplot(3,2,3); 
        lg = plot(squeeze(RF{n}(:,3,:)));% laser go
        set(lg, {'color'}, mat2cell(copper(nreps),ones(1,nreps)));
        ylim([min(min(squeeze(min(RF{n})))),max(max(squeeze(max(RF{n}))))]);
        slg.Title.String = 'lg';

        slng = subplot(3,2,4); 
        lng = plot(squeeze(RF{n}(:,4,:)));% laser nogo
        set(lng, {'color'}, mat2cell(copper(nreps),ones(1,nreps)));
        ylim([min(min(squeeze(min(RF{n})))),max(max(squeeze(max(RF{n}))))]);
        slng.Title.String = 'lng';
        
        
        sall = subplot(3,2,[5 6]); 
        hold on; plot(nanmean(RF{n}(:,1,:),3),'Color',[0 0.7 0 ]);% vis go
        hold on; plot(nanmean(RF{n}(:,2,:),3),'Color',[0.7 0 0 ]);% vis nogo
        hold on; plot(nanmean(RF{n}(:,3,:),3),'g');% laser go
        hold on; plot(nanmean(RF{n}(:,4,:),3),'r');% laser nogo
        legend({'vg','vn','lg','ln'})
        sall.Title.String = 'averaged over reps';
        
    end
end

function RF = calcERFforWindow(res,RF,nlags,ninputs,rep)
for n = 1:size(RF,2) % n is number of neurons
    for in = 1:ninputs
        for i=1:nlags
            RF{n}(i,in,rep) = res.params.model.C(n,:)*((res.params.model.A)^(i-1))*res.params.model.B(:,in);
        end
    end
end
end

% figure;
% h = plot(params.model.B);
% set(h, {'color'}, {[0 0.7 0 ]; [0.7 0 0 ]; [0 1 0]; [1 0 0]});
% legend({'vg','vn','lg','ln'})


