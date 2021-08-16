% loop over all animals, make plots and append

cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/matlabCode/scripts')


animallist ={'VL61','VL63','VL55','VL59',...
    'MPV33','MPV31','MPV34_2',...
    'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};


exptype = {'FF','FF','FF','FF',...
    'FF','FF','FF',...
    'FB','FB',...
    'FB','FB','FB','FB','FB'};

rootdir  = '/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/';
savedest = '/mnt/data/Mitra/cache/repos/figures/16082021/';
doSave = 1;
area = 'LM';

for animali = 1:length(animallist)
    animalname = animallist{animali};
    cd(fullfile(rootdir,animalname,'/trial_based/50msBins/'))
    PLDSresFiles = dir([area,'*.mat']);
    if length(PLDSresFiles) > 0
        % by default takes the latest file
        if length(PLDSresFiles)>1
            PLDSresFiles = PLDSresFiles(1);
        end
        res = load(PLDSresFiles.name);
        stimOn = sum(res.seq(1).u(1:2,:),1); % based on trial 1
        
        f = figure;f.Position = [159 519 1490 422];
        s1=subplot(1,3,1);
        for state=1:9
            hold on;plot(mean(cell2mat(arrayfun(@(x) x.posterior.xsm(state,:), res.seq,'UniformOutput',0)'),1),'b')
        end
        hold on;patch([find(stimOn),fliplr(find(stimOn))],[s1.YLim(1)+zeros(size(find(stimOn))),s1.YLim(2)+zeros(size(find(stimOn)))],...
            'y','FaceAlpha',0.1,'EdgeAlpha',0);
        s1.Title.String = 'state variables, averaged across all trials';
        s1.Title.FontWeight='normal';s1.Title.FontSize=11;
        xlabel('timebins');
        
        
        s2 = subplot(1,3,2);
        plot(res.seq(1).posterior.xsm','r') % single trial
        hold on;patch([find(stimOn),fliplr(find(stimOn))],[s2.YLim(1)+zeros(size(find(stimOn))),s2.YLim(2)+zeros(size(find(stimOn)))],...
            'y','FaceAlpha',0.1,'EdgeAlpha',0);
        s2.Title.String = 'states variables in an example single trial';
        s2.Title.FontWeight='normal';s2.Title.FontSize=11;
        xlabel('timebins');
        
        
        s3=subplot(1,3,3);
        for neuron = 1:size(res.seq(1).y,1)
            hold on;plot(mean(cell2mat(arrayfun(@(x) x.y(neuron,:), res.seq,'UniformOutput',0)'),1),'k')
        end
        hold on;patch([find(stimOn),fliplr(find(stimOn))],[s3.YLim(1)+zeros(size(find(stimOn))),s3.YLim(2)+zeros(size(find(stimOn)))],...
            'y','FaceAlpha',0.1,'EdgeAlpha',0);
        s3.Title.String = sprintf('spiking activity, averaged across all trials\n each trace is one neuron');
        s3.Title.FontWeight='normal';s3.Title.FontSize=11;
        xlabel('timebins');
        
        set(f.Children,'box','off')
        
        if doSave
            f.Name=[area,'-',animalname,'-50msBins']
            savename = fullfile(savedest,[area,'-',animalname,'-50msBins']);
            saveas(f,savename,'fig')
            cd(savedest)
        end
        
    end
    
    
end