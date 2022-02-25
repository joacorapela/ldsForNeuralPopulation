  figure; 
for r = 0:10
Mfile = dir(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/','VL66',...
             '//Joint_trial_based_splitContext_CntrlOnly/17msBins/*','RND',num2str(r),'*.mat']); % 17ms BaselineU = 0
cd(Mfile(end).folder)
res = load(Mfile(end).name);
v=res.varBound(find(~isnan(res.varBound)));
s= subplot(3,11,r+1);plot(res.varBound);s.Title.String =num2str(v(end));
subplot(3,11,r+11+1);imagesc(res.params.model.A);caxis([-0.5,0.5])
subplot(3,11,r+22+1);imagesc(res.params.model.B);caxis([-1.,1.])

end



%%
X = [];
for i = 1:length(seq)
    X(:,:,end+1) = seq(i).posterior.xsm;
end

figure;plot(nanmean(X,3)');legend
%% visualize u and y alignment
i = 1
figure;plot(seq(i).u');
hold on;plot(sum(seq(i).y(end-15:end,:))); % lm neurons
%%



heldoutN = 27
pred_bs = nan(length(seq),seq(1).T); % trials*timebins
yOrig_bs = nan(length(seq),seq(1).T); % trials*timebins
pred_ls = nan(length(seq),seq(1).T); % trials*timebins
yOrig_ls = nan(length(seq),seq(1).T); % trials*timebins

for tr = 1:length(seq)
   
    % for trial tr:
    z = params.model.C(heldoutN,:) *seq(tr).posterior.xsm + ...
        params.model.d(heldoutN);
    %if sum(sum(seq(tr).u(3:end,:),1))>0
        if sum(sum(seq(tr).u(2,:),1)) >0
           % disp('bla')
    pred_ls(tr,:) = exp(z);
    yOrig_ls(tr,:) = seq(tr).y(heldoutN,:);
       % end
    else
    pred_bs(tr,:) = exp(z);
    yOrig_bs(tr,:) = seq(tr).y(heldoutN,:);
    end
end


figure;


subplot(1,2,1);plot(nanmean(pred_bs,1));hold on;plot(nanmean(yOrig_bs,1));legend({'bs prediction','bs original trace'})
subplot(1,2,2);plot(nanmean(pred_ls,1));hold on;plot(nanmean(yOrig_ls,1));legend({'ls prediction','ls original trace'})
%%

g = load('Joint_PLDSfitRes_22_02_16_17_23_28_RND0_onlyCorrect_exGrooming_go.mat');
ng = load('Joint_PLDSfitRes_22_02_16_18_45_42_RND0_onlyCorrect_exGrooming_nogo.mat');
%%
animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};

for animali  = 1:6
    cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
    
    allA = nan(16,16,2);
    
    d = dir(['*RND0_onlyCorrect_exGrooming_go.mat']);
    load(d.name)
    allA(:,:,1) = params.model.A;
    d = dir(['*RND0_onlyCorrect_exGrooming_nogo.mat']);
    load(d.name)
    allA(:,:,2) = params.model.A;
    
    
    figure
    for lag = 1:2
        % LM pca
        
        [c,s,l]=pca(allA(1:8,9:16,lag));
        % not sure if transpose neede or not
        %hold on;plot(l); % almost 1d
        LMdir = c(:,1);
        
        B = allA(1:8,9:16,lag)*LMdir; % I think
        [u,v] = eig(allA(1:8,1:8,lag));
        [~,Ord]=sort(abs((diag(v))));
        proj = abs(inv(u) * B);
        hold on;plot(proj(Ord)./norm(proj))
    end
    legend
end
% how fast are V1 local dynamics in go vs nogo 
[u,v]= eig(g.params.model.A(1:8,1:8));
abs(diag(v)) 
imag(diag(v))

[u,v]= eig(ng.params.model.A(1:8,1:8));
abs(diag(v)) 
imag(diag(v))

% look at interareal matrix
figure;subplot(1,2,1);imagesc(g.params.model.A(1:8,9:16));caxis([-0.2,0.2])
subplot(1,2,2);imagesc(ng.params.model.A(1:8,9:16));caxis([-0.2,0.2])
 
% how low dim are modes in V1/LM?
[c,s,l]=pca(g.params.model.A(1:8,9:16)');l
[c,s,l]=pca(ng.params.model.A(1:8,9:16)');l
[c,s,l]=pca(g.params.model.A(1:8,9:16));l
[c,s,l]=pca(ng.params.model.A(1:8,9:16));l
%%
animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};

for animali  = 1:6
    cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
    
    allA = nan(16,16,2);
    
    d = dir(['*RND0_onlyCorrect_exGrooming_go.mat']);
    load(d.name)
    allA(:,:,1) = params.model.A;
    d = dir(['*RND0_onlyCorrect_exGrooming_nogo.mat']);
    load(d.name)
    allA(:,:,2) = params.model.A;
    
    
    figure
    for lag = 1:2
        % LM pca
        
        [c,s,l]=pca((allA(9:16,9:16,lag))');
        % not sure if transpose neede or not
        hold on;plot(l/norm(l)); % almost 1d
        LMdir = c(:,1);
        
        B = allA(1:8,9:16,lag)*LMdir; % I think
        [u,v] = eig(allA(1:8,1:8,lag));
        [~,Ord]=sort(abs((diag(v))));
        proj = abs(inv(u) * B);
       % hold on;plot(proj(Ord)./norm(proj))
    end
    legend
end




