% loop over cluster output 
% specs to print: area, animalname, splitdelays, bin size,inference method

% remove old for examining new results
cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/cluster/old/results')
animallist ={'VL61','VL63','VL55','VL59',...
    'MPV33','MPV31','MPV34_2',...
    'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66','MPV35_2'};

animali = 1;
binSizems = 50;
area = 'V1';


animalname = animallist{animali};
cd(sprintf('%s/trial_based_LONO_MS/%smsBins',animalname,num2str(binSizems)))

% by default reads tghe latest version
ll = nan(1,50);
for st = 1:numel(ll)
    File = dir(['nStates',num2str(st),'-','*.mat']);
    % if many, last one
    if numel(File) > 1
       File = File(end);
    end
    if numel(File)
        res = load(File.name);
        res = res.ModelSelection;
    end
end

