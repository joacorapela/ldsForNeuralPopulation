% CV (coefficient of variation) of laser and visual stimulus input receptive fields
% separated based on FF/FB and V1/LM

cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/matlabCode/scripts')

% animallist = {'VL61','VL63','VL55','VL59',...
%     'MPV33','MPV31','MPV34_2'};%,...
%    'MPV17','MPV18_2',...
%     'VL53','VL52','VL51','VL66'};%,'MPV35_2'};
animallist = { 'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};%,'MPV35_2'};
animalcolors = lines(7);

% set params
exptype = 'FB'; % set based on the selected animals
area = 'V1';
drawlines = 0;
binsize = '50ms/';  % 100
windowsize = '180'; % ms %300
AvWins = 1; % if 0 pools different fits of the same animal, if 1, averages.



%%% all 4 receptive fields (at zero delay) for laser anf
%%% visual stumuls input

Allvg = [];
Allvn = [];
Alllg = [];
Allln = [];
x = [];
y1 = [];
y2 = [];
for animali = 1%:length(animallist)
    animalname = animallist{animali};
    cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animalname,'/',binsize,windowsize])
    
    PLDSresFiles = dir([area,'*.mat']);
    nreps = length(PLDSresFiles);
    
    Animalvg = []; 
    Animalvn = []; 
    Animallg = []; 
    Animalln = []; 
    animalV = [];
    clear res
    for rep = 1:nreps
        res = load(PLDSresFiles(rep).name);
        vg = (res.params.model.C)*res.params.model.B(:,1); % eRF at delay 0 for all animals
        vn = (res.params.model.C)*res.params.model.B(:,2);
        lg = (res.params.model.C)*res.params.model.B(:,3);
        ln = (res.params.model.C)*res.params.model.B(:,4);
        
        Animalvg = [Animalvg,vg];
        Animalvn = [Animalvn,vn];
        Animallg = [Animallg,lg];
        Animalln = [Animalln,ln]; 
        animalV(end+1) = res.seq.posterior.varBound;
    end
   animalV = nanmean(animalV);
   x(end+1)=animalV;
    
%     Animalvg = mean(Animalvg,2);
%     Animalvn = mean(Animalvn,2);
%     Animallg = mean(Animallg,2);
%     Animalln = mean(Animalln,2);
    y1(end+1) = nanmean(std(Animalvg')'./abs(mean(Animalvg')')); % coef of variation
    y2(end+1) = nanmean(std(Animallg')'./abs(mean(Animallg')')); % coef of variation
    
    Animalvg = std(Animalvg')'./abs(mean(Animalvg')');
    Animalvn = std(Animalvn')'./abs(mean(Animalvn')');
    Animallg = std(Animallg')'./abs(mean(Animallg')');
    Animalln = std(Animalln')'./abs(mean(Animalln')');

    Allvg = [Allvg;(Animalvg)];
    Allvn = [Allvn;(Animalvn)];
    Alllg = [Alllg;(Animallg)];
    Allln = [Allln;(Animalln)];

end
figure;
s1 = subplot(1,1,1);
hold on;scatter(1+rand(1,length(Allvg))/3-0.5/3,Allvg,'.g');
hold on;scatter(2+rand(1,length(Allvn))/3-0.5/3,Allvn,'.r');
hold on;scatter(3+rand(1,length(Alllg))/3-0.5/3,Alllg,'.g');
hold on;scatter(4+rand(1,length(Allln))/3-0.5/3,Allln,'.r');

hold on;scatter(1,nanmedian(Allvg),'k+')
hold on;scatter(2,nanmedian(Allvn),'k+')
hold on;scatter(3,nanmedian(Alllg),'k+')
hold on;scatter(4,nanmedian(Allln),'k+')
xlim([0 5])
ylim([0 10])
s1.Title.String = [binsize,exptype,area,'-CV-RFs(delay0)'];
%%
% s2 = subplot(1,3,2);
% scatter(x,y)