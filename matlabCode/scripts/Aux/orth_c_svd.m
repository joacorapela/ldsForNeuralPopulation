function [CO,W] = orth_c_svd(params)

[U,S,V] = svd(params.model.C);
M = U*S;
N = V';

nV1= numel(find(params.model.CMask(:,1)));
nLM = numel(find(params.model.CMask(:,end)));
if nV1+nLM ~= size(params.model.C,1)
    error
end

ord = [find(sum(M(nV1+1:end,:),1) == 0),find(sum(M(nV1+1:end,:),1) ~= 0)];

%sum(sum(((M*N) - (M(:,ord)*N(ord,:))).^2))

%figure;imagesc(params.model.C - M*N)
%figure;imagesc(params.model.C - M(:,ord)*N(ord,:));

CO = M(:,ord);
W= N(ord,:);

%figure;imagesc(CO);
%cond(W)



end