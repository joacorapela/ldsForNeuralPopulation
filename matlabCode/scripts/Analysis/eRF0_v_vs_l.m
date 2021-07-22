% eRF_v (0) against eRF_l (0): comparing instantanious eRFs of visual vs
% laser inputs. (only 1 rep, no avg - to avoid removing relations by 
% averaging across fits)


cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/matlabCode/scripts')

animallist = {'VL61','VL63','VL55','VL59',...
    'MPV33','MPV31','MPV34_2'};%,...
% %    'MPV17','MPV18_2',...
% %    'VL53','VL52','VL51','VL66','MPV35_2'};
% animallist = { 'MPV17','MPV18_2',...
%     'VL53','VL52','VL51','VL66'};%,'MPV35_2'};
animalcolors = lines(7);

% set params
exptype = 'FF'; % set based on the selected animals
area = 'LM';
drawlines = 0;
binsize = '100ms/';  % 100
windowsize = '300'; % ms %300
AvWins = 1; % if 0 pools different fits of the same animal, if 1, averages.



%%% sign of the receptive fields of the cells (at zero delay) for laser vs.
%%% visual stumuls input
figure
Allgo_sign = [];
Allnogo_sign = [];
for animali = 1:length(animallist)
    animalname = animallist{animali};
    cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animalname,'/',binsize,windowsize])
    
    PLDSresFiles = dir([area,'*.mat']);
    nreps = length(PLDSresFiles);
    
    vg = nan(70,nreps);vn = nan(70,nreps);
    lg = nan(70,nreps);ln = nan(70,nreps);
    for rep = 6%1:nreps
        res = load(PLDSresFiles(rep).name);
        vg(1:size(res.params.model.C,1),rep) = (res.params.model.C)*res.params.model.B(:,1); % eRF at delay 0 for all animals
        vn(1:size(res.params.model.C,1),rep) = (res.params.model.C)*res.params.model.B(:,2);
        lg(1:size(res.params.model.C,1),rep) = (res.params.model.C)*res.params.model.B(:,3);
        ln(1:size(res.params.model.C,1),rep) = (res.params.model.C)*res.params.model.B(:,4);
    end
   % hold on; scatter(nanmean(vn,2),nanmean(ln,2),100,'k.')
    lns = nanmean(lg,2);
    vns = nanmean(vg,2);
    hold on; scatter(vns(find(lns<0)),lns(find(lns<0)),100,'b.')
    hold on; scatter(vns(find(lns>0)),lns(find(lns>0)),100,'r.')
    xlim([-1 3]); ylim([-0.6 1])
    xlabel('vn');ylabel('ln')
end

