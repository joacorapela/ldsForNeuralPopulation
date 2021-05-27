% similarity of laser and visual stimulus input weights (or receptive fields)
% separated based on FF/FB and V1/LM

cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/matlabCode/scripts')

animallist = {'VL61','VL63','VL55','VL59',...
    'MPV33','MPV31','MPV34_2'};%,...
%    'MPV17','MPV18_2',...
%    'VL53','VL52','VL51','VL66','MPV35_2'};
animalcolors = lines(7);

% set params
exptype = 'FF'; % set based on the selected animals
area = 'LM';
drawlines = 0;
binsize = '100ms/';
windowsize = '300'; % ms
AvWins = 1; % if 0 pools different fits of the same animal, if 1, averages.



%%% Angle of the weight vectors
Allgo_ang = [];
Allnogo_ang = [];
figure;s1 = subplot(1,3,1);
for animali = 1:length(animallist)
    animalname = animallist{animali};
    cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animalname,'/',binsize,windowsize])
    
    PLDSresFiles = dir([area,'*.mat']);
    nreps = length(PLDSresFiles);
    Animalgo = []; Animalnogo = [];
    for rep = 1:nreps
        res = load(PLDSresFiles(rep).name);
        Animalgo(end+1) = (res.params.model.B(:,1)./norm(res.params.model.B(:,1)))' * (res.params.model.B(:,3)./norm(res.params.model.B(:,3)));
        Animalnogo(end+1) = (res.params.model.B(:,2)./norm(res.params.model.B(:,2)))' * (res.params.model.B(:,4)./norm(res.params.model.B(:,4)));
    end
    if AvWins
        Animalgo = mean(Animalgo);
        Animalnogo = mean(Animalnogo);
    end
    hold on;scatter(ones(1,length(Animalgo)),Animalgo,'g');
    hold on;scatter(1+ones(1,length(Animalnogo)),Animalnogo,'r');
    if drawlines
        for i=1:length(Animalgo)
            hold on;line([1,2],[Animalgo(i),Animalnogo(i)],'Color',animalcolors(animali,:))
        end
    end
    xlim([0 3])
    Allgo_ang = [Allgo_ang Animalgo];
    Allnogo_ang = [Allnogo_ang Animalnogo];
end
hold on;scatter(1,nanmean(Allgo_ang),'k+')
hold on;scatter(2,nanmean(Allnogo_ang),'k+')
s1.Title.String = [exptype,area,'-Ang'];

%%% sign of the receptive fields of the cells (at zero delay) for laser vs.
%%% visual stumuls input
s2 = subplot(1,3,2);
Allgo_sign = [];
Allnogo_sign = [];
for animali = 1:length(animallist)
    animalname = animallist{animali};
    cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animalname,'/',binsize,windowsize])
    
    PLDSresFiles = dir([area,'*.mat']);
    nreps = length(PLDSresFiles);
    
    Animalgo_sign = []; Animalnogo_sign = [];
    for rep = 1:nreps
        res = load(PLDSresFiles(rep).name);
        vg = (res.params.model.C)*res.params.model.B(:,1); % eRF at delay 0 for all animals
        vn = (res.params.model.C)*res.params.model.B(:,2);
        lg = (res.params.model.C)*res.params.model.B(:,3);
        ln = (res.params.model.C)*res.params.model.B(:,4);
        
%        should they be normalized?
%        Animalgo_sign = [Animalgo_sign;(vg/norm(vg)).* (lg/norm(lg))];
%        Animalnogo_sign = [Animalnogo_sign;(vn/norm(vn)).* (ln/norm(ln))];
       
       Animalgo_sign = [Animalgo_sign;vg.*lg];
       Animalnogo_sign = [Animalnogo_sign;vn.*ln];
       
%        try this:
%        Animalgo_sign = [Animalgo_sign;vg-vn];
%        Animalnogo_sign = [Animalnogo_sign;lg-ln];
       
    end
    if AvWins
        Animalgo_sign = mean(Animalgo_sign);
        Animalnogo_sign = mean(Animalnogo_sign);
    end
    hold on;scatter(ones(1,length(Animalgo_sign)),Animalgo_sign,'g');
    hold on;scatter(1+ones(1,length(Animalnogo_sign)),Animalnogo_sign,'r');
    xlim([0 3]);
    Allgo_sign = [Allgo_sign;Animalgo_sign];
    Allnogo_sign = [Allnogo_sign;Animalnogo_sign];
end
hold on;scatter(1,nanmean(Allgo_sign),'k+')
hold on;scatter(2,nanmean(Allnogo_sign),'k+')
s2.Title.String = [exptype,area,'-RFsign'];

%%% all 4 receptive fields (at zero delay) for laser anf
%%% visual stumuls input
s3 = subplot(1,3,3);
Allvg = [];
Allvn = [];
Alllg = [];
Allln = [];
for animali = 1:length(animallist)
    animalname = animallist{animali};
    cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animalname,'/',binsize,windowsize])
    
    PLDSresFiles = dir([area,'*.mat']);
    nreps = length(PLDSresFiles);
    
    Animalvg = []; 
    Animalvn = []; 
    Animallg = []; 
    Animalln = []; 
    for rep = 1:nreps
        res = load(PLDSresFiles(rep).name);
        vg = (res.params.model.C)*res.params.model.B(:,1); % eRF at delay 0 for all animals
        vn = (res.params.model.C)*res.params.model.B(:,2);
        lg = (res.params.model.C)*res.params.model.B(:,3);
        ln = (res.params.model.C)*res.params.model.B(:,4);
        
        Animalvg = [Animalvg;vg];
        Animalvn = [Animalvn;vn];
        Animallg = [Animallg;lg];
        Animalln = [Animalln;ln]; 
    end
    if AvWins
        Animalvg = mean(Animalvg);
        Animalvn = mean(Animalvn);
        Animallg = mean(Animallg);
        Animalln = mean(Animalln);
    end
    Allvg = [Allvg;Animalvg];
    Allvn = [Allvn;Animalvn];
    Alllg = [Alllg;Animallg];
    Allln = [Allln;Animalln];
end
hold on;scatter(ones(1,length(Allvg)),Allvg,'g');
hold on;scatter(1+ones(1,length(Allvn)),Allvn,'r');
hold on;scatter(2+ones(1,length(Alllg)),Alllg,'g');
hold on;scatter(3+ones(1,length(Allln)),Allln,'r');

hold on;scatter(1,nanmean(Allvg),'k+')
hold on;scatter(2,nanmean(Allvn),'k+')
hold on;scatter(3,nanmean(Alllg),'k+')
hold on;scatter(4,nanmean(Allln),'k+')
xlim([0 5])
s3.Title.String = [exptype,area,'-RFs(delay0)'];