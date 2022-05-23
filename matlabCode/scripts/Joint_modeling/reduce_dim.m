function TrYcntrl_red = reduce_dim(TrYcntrl,ndims,showplot)

X = reshape(TrYcntrl,[],size(TrYcntrl,3));%timepointsandtrials*neurons

[coefs,scores,vars] = pca(X,'Centered',0); %% [TODO] check if it should be centered or not
% X = scores*coefs'

[~,order] = sort(vars,'descend');
coefs = coefs(:,order);

coefs_red = [coefs(:,1:ndims),zeros(size(coefs,1),size(coefs,2)-ndims)]; % mean subtraction?
Xred = scores*coefs_red';

%figure;imagesc(X);figure;imagesc(Xred)
TrYcntrl_red=reshape(Xred,size(TrYcntrl,1),size(TrYcntrl,2),size(TrYcntrl,3));
if showplot
    figure;imagesc(squeeze(nanmean(TrYcntrl(:,:,:),1)));figure;imagesc(squeeze(nanmean(TrYcntrl_red(:,:,:),1)))
    figure;plot(nanmean(TrYcntrl(:,:,23),1));hold on;plot(nanmean(TrYcntrl_red(:,:,23),1),'r');legend({'original','low dim'})
    % figure;plot(TrYcntrl(2,:,6));hold on;plot(TrYcntrl_red(2,:,6),'r');
end


% %% PCA scratch
% %X = (squeeze(nanmean(TrY(:,timewindow,:),1)));%timepoints*neurons
% X = reshape(TrY,[],size(TrY,3));%timepointsandtrials*neurons
% 
% [coefs,scores,~] = pca(X,'Centered',0);
% % X = scores*coefs'
% coefs_red = [coefs(:,1:16),zeros(size(coefs,1),size(coefs,2)-16)]; % assuming descending order, mean subtraction?
% Xred = scores*coefs_red';
% %figure;imagesc(X);figure;imagesc(Xred)
% 
% TrYOut=reshape(Xred,size(TrY,1),size(TrY,2),size(TrY,3));
% figure;imagesc(squeeze(nanmean(TrY(:,:,:),1)));figure;imagesc(squeeze(nanmean(TrYOut(:,:,:),1)))
% 
% %%%%%
% X = reshape(TrYcntrl,[],size(TrY,3));%timepointsandtrials*neurons
% 
% [coefs,scores,~] = pca(X,'Centered',0);
% % X = scores*coefs'
% coefs_red = [coefs(:,1:16),zeros(size(coefs,1),size(coefs,2)-16)]; % assuming descending order, mean subtraction?
% Xred = scores*coefs_red';
% %figure;imagesc(X);figure;imagesc(Xred)
% 
% TrYcntrl_red=reshape(Xred,size(TrYcntrl,1),size(TrYcntrl,2),size(TrYcntrl,3));
% figure;imagesc(squeeze(nanmean(TrYcntrl(:,:,:),1)));figure;imagesc(squeeze(nanmean(TrYcntrl_red(:,:,:),1)))