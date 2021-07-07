% similarity of laser and visual stimulus input weights (or receptive fields)
% separated based on FF/FB and V1/LM

cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/matlabCode/scripts')

% animallist = {'VL61','VL63','VL55','VL59',...
%     'MPV33','MPV31','MPV34_2'};%,...
% %    'MPV17','MPV18_2',...
% %    'VL53','VL52','VL51','VL66','MPV35_2'};
animallist = { 'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};%,'MPV35_2'};
animalcolors = lines(7);

% set params
exptype = 'FB'; % set based on the selected animals
area = 'V1';
drawlines = 0;
binsize = '100ms/';  % 100
windowsize = '300'; % ms %300
AvWins = 1; % if 0 pools different fits of the same animal, if 1, averages.



%%% Angle of the weight vectors
Allgo_ang = [];
Allnogo_ang = [];
figure;%s1 = subplot(1,3,1);
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
%     hold on;scatter(1+rand(1,length(Animalgo))/5-0.1,Animalgo,'.g');
%     hold on;scatter(2+rand(1,length(Animalnogo))/5-0.1,Animalnogo,'.r');
    if drawlines
        for i=1:length(Animalgo)
            hold on;line([1,2],[Animalgo(i),Animalnogo(i)],'Color',animalcolors(animali,:))
        end
    end
    xlim([0 3])
    Allgo_ang = [Allgo_ang Animalgo];
    Allnogo_ang = [Allnogo_ang Animalnogo];
end
% hold on;scatter(1,nanmean(Allgo_ang),'k+')
% hold on;scatter(2,nanmean(Allnogo_ang),'k+')
% s1.Title.String = [binsize,exptype,area,'-Ang'];

%%% sign of the receptive fields of the cells (at zero delay) for laser vs.
%%% visual stumuls input
s2 = subplot(1,2,1);
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
       
       Animalgo_sign = [Animalgo_sign,vg.*lg]; % , ;
       Animalnogo_sign = [Animalnogo_sign,vn.*ln];
       
%        try this:
%        Animalgo_sign = [Animalgo_sign;vg-vn];
%        Animalnogo_sign = [Animalnogo_sign;lg-ln];
       
    end
    if AvWins
        Animalgo_sign = mean(Animalgo_sign,2);% 2 nan
        Animalnogo_sign = mean(Animalnogo_sign,2);
    else
        Animalgo_sign = reshape(Animalgo_sign,[],1);
        Animalnogo_sign = reshape(Animalnogo_sign,[],1);
    end
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


%%% all 4 receptive fields (at zero delay) for laser anf
%%% visual stumuls input
s3 = subplot(1,2,2);
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
    repV = [];
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
        repV(end+1) = res.seq.posterior.varBound;
    end
    if 0
       Animalvg = cleanupCI(Animalvg,nreps);
       Animalvn = cleanupCI(Animalvn,nreps);
       Animallg = cleanupCI(Animallg,nreps);
       Animalln = cleanupCI(Animalln,nreps);
    end
    if 0
        [~,inds] = max(repV);
        Animalvg = Animalvg(:,inds);
        Animalvn = Animalvn(:,inds);
        Animallg = Animallg(:,inds);
        Animalln = Animalln(:,inds);
    end
    if AvWins
        Animalvg = mean(Animalvg,2);
        Animalvn = mean(Animalvn,2);
        Animallg = mean(Animallg,2);
        Animalln = mean(Animalln,2);
    else
        Animalvg = reshape(Animalvg,[],1);
        Animalvn = reshape(Animalvn,[],1);
        Animallg = reshape(Animallg,[],1);
        Animalln = reshape(Animalln,[],1);
    end
    Allvg = [Allvg;Animalvg];
    Allvn = [Allvn;Animalvn];
    Alllg = [Alllg;Animallg];
    Allln = [Allln;Animalln];
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

s3.Title.String = [binsize,exptype,area,'-RFs(delay0)'];
s3.XTick = [1 2 3 4];
s3.XTickLabel = {'vg','vn','lg','ln'};

function Animalvg = cleanupCI(Animalvg,nrep)
 for rows = 1:size(Animalvg,1)
            if (nanmean(Animalvg(rows,:)))> 0 
                if (nanmean(Animalvg(rows,:)) - 2*nanstd(Animalvg(rows,:))/sqrt(nrep))<0
                    Animalvg(rows,:) = 0;
                end

            elseif (nanmean(Animalvg(rows,:)))< 0 
                if (nanmean(Animalvg(rows,:)) + 2*nanstd(Animalvg(rows,:))/sqrt(nrep))>0
                    Animalvg(rows,:) = 0;
                end
            end
 end
end