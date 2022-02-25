% after running scratch_V1NonlinearReadout
clc
numanimals=6;
numtimepoints = 8;
Mbs= nan(1,numtimepoints);
Mls= nan(1,numtimepoints);
Ebs= nan(1,numtimepoints);
Els= nan(1,numtimepoints);

allgocells = cell(1,numtimepoints);
allnogocells = cell(1,numtimepoints);


for timepoint = 1:numtimepoints
    
    allgocells{timepoint} = [];
    allnogocells{timepoint} = [];
    
    for animalnum=1:numanimals
              
        targetcell = V1cells(find(cellfun(@(x) x.simulcode,V1cells) == animalnum ));
        % ntrials*ncells:
        gobs = cell2mat(cellfun(@(x) x.laAbs.go{timepoint},targetcell,'UniformOutput',0));
        nogobs = cell2mat(cellfun(@(x) x.laAbs.nogo{timepoint},targetcell,'UniformOutput',0));
        
        bsall = [gobs;nogobs];
        
        
        % equalize trials
        gobs(:,isnan(sum(bsall,1))) = [];
        nogobs(:,isnan(sum(bsall,1))) = [];
        
        mintrialnum = min([size(gobs,1),size(nogobs,1)]);
        gobs = gobs(1:mintrialnum,:);
        nogobs = nogobs(1:mintrialnum,:);
        
        allgocells{timepoint} = [ allgocells{timepoint}, nanmean(gobs,1)];
        allnogocells{timepoint} = [ allnogocells{timepoint}, nanmean(nogobs,1)];
    end
    
    allgocells{timepoint}(end+1:100) = nan;
    allnogocells{timepoint}(end+1:100) = nan;
    
    Mgo(timepoint) =  nanmedian(allgocells{timepoint}); 
    Mnogo(timepoint) =  nanmedian(allnogocells{timepoint});  
    Mboth(timepoint) =  nanmedian([allgocells{timepoint},allnogocells{timepoint}]); 
 
    Ego(timepoint) =  nanstd(allgocells{timepoint}); 
    Enogo(timepoint) =  nanstd(allnogocells{timepoint});  
    Eboth(timepoint) =  nanstd([allgocells{timepoint},allnogocells{timepoint}]); 
    
end
%%
figure;boxplot(cell2mat(allgocells')'*1000/70,'PlotStyle','traditional','Symbol','','Colors','g',...
    'MedianStyle','line','Widths',0.5,'BoxStyle','outline') % **** Hz conversion only for an70
ylabel('firing rate (Hz)')
figure;boxplot(cell2mat(allnogocells')'*1000/70,'PlotStyle','traditional','Symbol','','Colors','r',...
    'MedianStyle','line','Widths',0.5,'BoxStyle','outline') % **** Hz conversion only for an70
ylabel('firing rate (Hz)')
figure;boxplot(0.5*(cell2mat(allgocells')'+cell2mat(allgocells')')*1000/70,'PlotStyle','traditional','Symbol','','Colors','k',...
    'MedianStyle','line','Widths',0.5,'BoxStyle','outline') % **** Hz conversion only for an70
ylabel('firing rate (Hz)')
%%
figure;shadedErrorBar([],Mgo,Ego,'g',1);
hold on;
shadedErrorBar([],Mnogo,Enogo,'r',1)

figure;shadedErrorBar([],Mboth,Eboth,'k',1);


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
