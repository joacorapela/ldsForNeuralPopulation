


figure;
allA = nan(10,10,8);

for i = 0:7
    
    d = dir(['JointSep_lag',num2str(i),'*.mat']);
    load(d.name)
    subplot(1,8,i+1)
    imagesc((params.model.A));
    caxis([-0.5, 0.5])
    %colorbar
    allA(:,:,i+1) = params.model.A;
end
    
%%
figure;cmap =jet(8);
for i = 1:8 % 5*5 A for V1 is normal only at  lag 2 and 5! and
 [u,v]=eig(allA(1:5,1:5,i));
 abs(v);
  hold on;plot(sort(abs(diag(v))),'Color',cmap(i,:))
end
legend

%% 
figure;cmap =jet(8);
for i = 1:8 
 [u,v]=eig(allA(1:5,6:10,i));
 hold on;plot(sort(abs(diag(v))),'Color',cmap(i,:))
 %pause(1)
end
legend
%%
figure;plot(squeeze(allA(3,4,:)))
hold on;plot(squeeze(allA(4,3,:)))
%% eigen decompose: magnitude of lambda on one axis and separation on the other
figure;
for i = 1:8
 [u,v]=eig(allA(:,:,i));
 
 V1norm = vecnorm(abs(u(1:5,:)));
 LMnorm = vecnorm(abs(u(6:10,:)));
 sep = ((V1norm-LMnorm)./(V1norm+LMnorm));
% subplot(1,8,i);histogram(sep,-1:0.2:1)
 subplot(1,8,i);scatter(sep,abs(diag(v)));xlim([-1,1]);ylim([0,1]);
 corrcoef(sep,abs(diag(v)))
end
%% eigen decompose and order eigen vecs

figure;
for i = 1:8
 [u,v]=eig(allA(:,:,i));
 [~,Ord] = sort(abs(diag(v)));
 subplot(1,8,i);imagesc(abs(u(:,Ord)))
end
%% project in LM
figure;
for i = 1:8
    subplot(1,8,i)
    for j = 6:10
    hold on;plot(allA(1:5,j,i))
    end
end
%% project in V1
figure;
for i = 1:8
    subplot(1,8,i)
    for j = 1:5
    hold on;plot(allA(j,6:10,i))
    end
end
%% 
% for example: this is the LM mode that influences V1 the most: figure;plot(allA(1:5,8,2))

%% how much intra mode is aligned with inter (is LM->V1 dyn is trated differently than just a B vec?)
% consistency across animals coud also e checked
% instead of each LM mode, it could also be eig of LM, 

m = 3;% LM mode number. We go through all lets say

figure
for lag = 1:8
B = allA(1:5,5+m,lag);
[u,v] = eig(allA(1:5,1:5,lag));
[~,Ord]=sort(abs(diag(v)));
proj = abs(inv(u) * B);
hold on;plot(proj(Ord)./norm(proj))
abs(diag(v))
end
legend


figure
for m = 1:5
for lag = 1:8
B = allA(1:5,5+m,lag);
[u,v] = eig(allA(1:5,1:5,lag));
[~,Ord]=sort(abs(diag(v)));
proj = abs(inv(u) * B);
subplot(1,8,lag); hold on;plot(proj(Ord)./norm(proj))
ylim([0 1])
end
end
legend




figure
for lag = 1:8
B = allA(1:5,5+3,lag);
[u,v] = eig(allA(1:5,1:5,lag));
[SortedV,Ord]=sort(abs(diag(v)));
proj = abs(inv(u) * B);
subplot(1,8,lag);scatter(SortedV,proj(Ord)./norm(proj))

end

%%%%%  eig of LM
%%*******************************************
% not really : not eig but some sort of pca! It seems LM influence on V1 is
% almost 1D - lets find this D use iy as LMdir


figure
for lag = 1:8
    % LM pca
    
    [c,s,l]=pca(allA(1:5,6:10,lag));
    % not sure if transpose neede or not
    %figure;plot(l); % almost 1d
    LMdir = c(:,1);
  
    B = allA(1:5,6:10,lag)*LMdir; % I think
    [u,v] = eig(allA(1:5,1:5,lag));
    [~,Ord]=sort(abs(diag(v)));
    proj = abs(inv(u) * B);
    hold on;plot(proj(Ord)./norm(proj))
end
legend

%% how much V1-V1 (local) is aligned with v1-lm-v1(one loop to lm and back)
figure;subplot(1,2,1);imagesc(allA(1:5,1:5,2));subplot(1,2,2);imagesc(allA(1:5,6:10,2));

for i = 1:8
figure;subplot(1,2,1);imagesc(allA(1:5,1:5,i));subplot(1,2,2);imagesc(allA(1:5,6:10,i)*allA(6:10,1:5,i));caxis([-0.5, 0.5])
[c,s,l]=pca((allA(1:5,6:10,i)*allA(6:10,1:5,i)));%figure;plot(l) %%%% This is literaly id without transpose. Does it need a transpose?
subspace(allA(1:5,1:5,i)',(allA(1:5,6:10,i)*allA(6:10,1:5,i))')
end
% similarity of the 2 vecs needs to be tracked  -- I dont think it is
% actually subspace

figure
for i = 1:8
[c,s,l]=pca((allA(1:5,6:10,i)*allA(6:10,1:5,i)));%figure;plot(l) %%%% This is literaly id without transpose. Does it need a transpose?

 B = c(:,1);
  [u,v] = eig(allA(1:5,1:5,i));
        [~,Ord]=sort(abs(diag(v)));
        proj = abs(inv(u) * B);
        hold on;plot(proj(Ord)./norm(proj))
end
legend
%% comparinf the dimensionality of intra and inter subspaces - private dimension stuff. Is it true for all animals?
figure;
%%% **** Is transpose needed or not? does it change over time?
for lag = 1:8
[c,s,l]=pca(allA(1:5,1:5,lag));
subplot(1,8,lag);hold on;plot(l/sum(l))

[c,s,l]=pca(allA(1:5,6:10,lag)); % try both this and transpose: How low dim is influence on V1 or from LM? (on V1 is kinda by definition 1d? or not?)
subplot(1,8,lag);hold on;plot(l/sum(l))

[c,s,l]=pca(allA(6:10,6:10,lag));
subplot(1,8,lag);hold on;plot(l/sum(l))

legend
end
%% comparing pca(A3) wih pca(A3') - prominent mode in V1 and LM?
figure
dirLM = nan(5,8);
dirV1 = nan(5,8);
pcn = 1;
for lag = 1:8
    % LM pca
    
    [c,s,l]=pca(allA(1:5,6:10,lag));
    % not sure if transpose neede or not
    %figure;plot(l); % almost 1d
    dirLM(:,lag) = c(:,pcn);
    
     [c,s,l]=pca(allA(1:5,6:10,lag)');
    % not sure if transpose neede or not
    %figure;plot(l); % almost 1d
    dirV1(:,lag) = c(:,pcn);
    
end
figure;subplot(2,2,1);imagesc(dirLM);subplot(2,2,2);imagesc(dirV1);
subplot(2,2,3);plot(dirLM);subplot(2,2,4);plot(dirV1);

%% more animals
animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66','MPV35_2'};

for animali  = 1:7
    cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_split_0/25msBins'])
    
    allA = nan(10,10,8);
    for i = 0:7
        
        d = dir(['JointSep_lag',num2str(i),'*.mat']);
        load(d.name)
        allA(:,:,i+1) = params.model.A;
    end
    
    figure
    for lag = 1:8
        % LM pca
        
        [c,s,l]=pca(allA(1:5,6:10,lag));
        % not sure if transpose neede or not
        %figure;plot(l); % almost 1d
        LMdir = c(:,1);
        
        B = allA(1:5,6:10,lag)*LMdir; % I think
        [u,v] = eig(allA(1:5,1:5,lag));
        [~,Ord]=sort(abs(diag(v)));
        proj = abs(inv(u) * B);
        hold on;plot(proj(Ord)./norm(proj))
    end
    legend
end
