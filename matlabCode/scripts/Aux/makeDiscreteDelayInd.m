function array2 = makeDiscreteDelayInd(array)
% returns an array with the same size as input (LaserDelayBinned). Its
% elements can take the following values: nan for no laser trials, 0-7 for
% laserdelays 0 to 7, and 8 for laserdelays that don't fit into any
% discrete delay values, these trials shoulld be discarded
%array = res.LaserDelayBinned;

freqs = nan(1,max(array));
for i = 1:max(array)
    freqs(i) = numel(find(array==i));
end

[val,ind] = sort(freqs,'descend');

Lags = ind(1:8);
Lags = sort(Lags);
Discard = ind(9:end);

array2 = array;
for i = 1:length(Lags)
    array2(array == Lags(i)) = i-1;
end
for i = 1:length(Discard)
    array2(array == Discard(i)) = 8;
end