
% alignment of lambdas (eigen values of dynamics matrix) and alphas
% (projection of input weights onto the eigen vectors of tyhe dynamics
% matrix)

animalname = 'MPV17';% VL53 fb %VL59 ff
binsize= '50ms'; % or '' for old files
area = 'V1';

real_D = [];
alpha1 = [];
alpha2 = [];
alpha3 = [];
alpha4 = [];


for rep = 1 % try also not in order but overlayed: order alphas then plot them based on index on real_D

trainDurSecs = 180;%180
analysisStartTimeSecs = (rep)*trainDurSecs; % choose the model to load here (start time)
cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results')
cd(fullfile(animalname,binsize))

% load timeseries file
% load PLDS results
cd(num2str(trainDurSecs))
PLDSresFiles = dir([area,sprintf('*_%s.mat',num2str(analysisStartTimeSecs))]);
res = load(PLDSresFiles.name);

[V,D] = eig(res.params.model.A);
% do real and abs
real_D = [real_D, abs(diag(D))'];% chnge ' 
alpha= nan(4,size(V,2));
for n = 1:size(V,2)
    alpha(1,n) = res.params.model.B(:,1)'*V(:,n);
    alpha(2,n) = res.params.model.B(:,2)'*V(:,n);
    alpha(3,n) = res.params.model.B(:,3)'*V(:,n);
    alpha(4,n) = res.params.model.B(:,4)'*V(:,n);
end
norm = 0
if norm
    alpha1 = [alpha1, (abs(alpha(1,:)) - min(abs(alpha(1,:))))/(max(abs(alpha(1,:))) - min(abs(alpha(1,:))))];
    alpha2 = [alpha2,(abs(alpha(2,:)) - min(abs(alpha(2,:))))/(max(abs(alpha(2,:))) - min(abs(alpha(2,:))))];
    alpha3 = [alpha3, (abs(alpha(3,:)) - min(abs(alpha(3,:))))/(max(abs(alpha(3,:))) - min(abs(alpha(3,:))))];
    alpha4 = [alpha4, (abs(alpha(4,:)) - min(abs(alpha(4,:))))/(max(abs(alpha(4,:))) - min(abs(alpha(4,:))))];
else
    alpha1 = [alpha1, alpha(1,:)];
    alpha2 = [alpha2, alpha(2,:)];
    alpha3 = [alpha3, alpha(3,:)];
    alpha4 = [alpha4, alpha(4,:)];   
end
end

[~,ord]=sort(real_D);
figure;plot(real_D(ord),'k');
hold on;plot(1*abs(alpha1(ord)),'Color',[0 0.6 0])
hold on;plot(1*abs(alpha2(ord)),'Color',[0.6 0 0])
hold on;plot(1*abs(alpha3(ord)),'Color',[0 1 0])
hold on;plot(1*abs(alpha4(ord)),'Color',[1 0 0]);

figure;
hold on;scatter(real_D(ord),abs(alpha1(ord)),'+','MarkerEdgeColor',[0 0.6 0])
hold on;scatter(real_D(ord),abs(alpha2(ord)),'+','MarkerEdgeColor',[0.6 0 0])
hold on;scatter(real_D(ord),abs(alpha3(ord)),'+','MarkerEdgeColor',[0 1 0])
hold on;scatter(real_D(ord),abs(alpha4(ord)),'+','MarkerEdgeColor',[1 0 0])

% hold on;plot(10*real(alpha(1,ord)),'Color',[0 0.6 0])
% hold on;plot(10*real(alpha(2,ord)),'Color',[0.6 0 0])
% hold on;plot(10*real(alpha(3,ord)),'Color',[0 1 0])
% hold on;plot(10*real(alpha(4,ord)),'Color',[1 0 0])


% use diag(d) as color
%%%

% cmap = jet(9);
% figure
% 
% for n=1:9
%     real_a = real(alpha(1,n));
%     imag_a = imag(alpha(1,n));
%     
%     hold on;plot(real_a,imag_a,'+','Color',cmap(ord(n),:))
% end
