% similarity of laser and visual stimulus input weights (or receptive fields)
% separated based on FF/FB and V1/LM

cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/matlabCode/scripts')
addpath('/mnt/data/Mitra/cache/codes/export_fig')
% set params
exptype = 'FB'; % set based on the selected animals
area = 'LM';
binsize = '50';  % 100
onefig = 1;
setylim = 1; % truncates y axis at predefined lims. use only when range is known

if strcmp(exptype,'FF')
    animallist = {'VL61','VL63','VL55','VL59',...
        'MPV33','MPV31','MPV34_2'};%,...
elseif strcmp(exptype,'FB')
    animallist = { 'MPV17','MPV18_2',...
        'VL53','VL52','VL51','VL66'};%,'MPV35_2'};
end

figure;
if ~onefig
    %%% Angle of the weight vectors
    Allgo_ang = [];
    Allnogo_ang = [];
    for animali = 1:length(animallist)
        animalname = animallist{animali};
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animalname,...
            '/trial_based/',binsize,'msBins/'])
        
        PLDSresFiles = dir([area,'*.mat']);
        if length(PLDSresFiles)>1
            PLDSresFiles = PLDSresFiles(1);
        end
        res = load(PLDSresFiles.name);
        
        Animalgo(end+1) = (res.params.model.B(:,1)./norm(res.params.model.B(:,1)))' * (res.params.model.B(:,3)./norm(res.params.model.B(:,3)));
        Animalnogo(end+1) = (res.params.model.B(:,2)./norm(res.params.model.B(:,2)))' * (res.params.model.B(:,4)./norm(res.params.model.B(:,4)));
        
        
        xlim([0 3])
        Allgo_ang = [Allgo_ang Animalgo];
        Allnogo_ang = [Allnogo_ang Animalnogo];
    end
    %%% sign of the receptive fields of the cells (at zero delay) for laser vs.
    %%% visual stumuls input
    s2 = subplot(1,2,1);
    Allgo_sign = [];
    Allnogo_sign = [];
    for animali = 1:length(animallist)
        animalname = animallist{animali};
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animalname,...
            '/trial_based/',binsize,'msBins/'])
        PLDSresFiles = dir([area,'*.mat']);
        if length(PLDSresFiles)>1
            PLDSresFiles = PLDSresFiles(1);
        end
        res = load(PLDSresFiles.name);
        vg = (res.params.model.C)*res.params.model.B(:,1); % eRF at delay 0 for all animals
        vn = (res.params.model.C)*res.params.model.B(:,2);
        lg = (res.params.model.C)*res.params.model.B(:,3);
        ln = (res.params.model.C)*res.params.model.B(:,4);
        Animalgo_sign = vg.*lg; % , ;
        Animalnogo_sign = vn.*ln;
        
        hold on;scatter(1+rand(1,length(Animalgo_sign))/5-0.1,Animalgo_sign,'.g');
        hold on;scatter(2+rand(1,length(Animalnogo_sign))/5-0.1,Animalnogo_sign,'.r');
        xlim([0 3]);
        ylim([-2,1]);%remove
        Allgo_sign = [Allgo_sign;Animalgo_sign];
        Allnogo_sign = [Allnogo_sign;Animalnogo_sign];
    end
    hold on;scatter(1,nanmean(Allgo_sign),'k+')
    hold on;scatter(2,nanmean(Allnogo_sign),'k+')
    s2.Title.String = [binsize,exptype,area,'-RFsign'];
    s2.XTick = [1 2];
    s2.XTickLabel = {'vg.lg','vn.ln'};
    s3 = subplot(1,1,1);
else
    s3 = subplot(1,1,1);
end
%%% all 4 receptive fields (at zero delay) for laser anf
%%% visual stumuls input

Allvg = [];
Allvn = [];
Alllg = [];
Allln = [];
for animali = 1:length(animallist)
    animalname = animallist{animali};
    cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animalname,...
        '/trial_based/',binsize,'msBins/'])
    PLDSresFiles = dir([area,'*.mat']);
    if length(PLDSresFiles)>1
        PLDSresFiles = PLDSresFiles(1);
    end
    try
        res = load(PLDSresFiles.name);
        vg = (res.params.model.C)*res.params.model.B(:,1); % eRF at delay 0 for all animals
        vn = (res.params.model.C)*res.params.model.B(:,2);
        lg = (res.params.model.C)*res.params.model.B(:,3);
        ln = (res.params.model.C)*res.params.model.B(:,4);
        
        Animalvg = vg;
        Animalvn = vn;
        Animallg = lg;
        Animalln = ln;
        
        Allvg = [Allvg;Animalvg];
        Allvn = [Allvn;Animalvn];
        Alllg = [Alllg;Animallg];
        Allln = [Allln;Animalln];
    catch
        disp('skipped animal')
    end
end
hold on;scatter(1+rand(1,length(Allvg))/5-0.1,Allvg,'g.');
hold on;scatter(2+rand(1,length(Allvn))/5-0.1,Allvn,'r.');
hold on;scatter(3+rand(1,length(Alllg))/5-0.1,Alllg,'g.');
hold on;scatter(4+rand(1,length(Allln))/5-0.1,Allln,'r.');

hold on;scatter(1,nanmean(Allvg),'k+')
hold on;scatter(2,nanmean(Allvn),'k+')
hold on;scatter(3,nanmean(Alllg),'k+')
hold on;scatter(4,nanmean(Allln),'k+')
xlim([0 5])
if setylim
    ylim([-4,4])
end
hold on; line([0,5],[0,0],'Color','k','LineStyle','--')

s3.Title.String = ['eRFs(0)-',exptype,'-',area,'-',binsize,'ms'];
s3.XTick = [1 2 3 4];
s3.XTickLabel = {'Vg','Vn','Lg','Ln'};
s3.FontSize = 14;

f=gcf;
set(f,'Color','w');