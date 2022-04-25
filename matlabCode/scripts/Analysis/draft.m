
[~,m]=sort(sum(abs(params.model.A(1:5,6:10)),1),'descend')

%figure;
x = nan(5,40);
for i = 1:5
%subplot(length(m),1,i);
x(i,:)=(mean(cell2mat(arrayfun(@(x) x.posterior.xsm(5+i,:),seq,'UniformOutput',0)'),1));
end

figure
for i=1:5
    subplot(5,1,i);
    hold on; plot(params.model.A(1:5,5+i)'*x)
end

figure
for i=1:5
    subplot(5,1,i);
    hold on;plot(x(i,:))
end
%%
norms = nan(1,18);
for i = 1:18
norms(i)=norm(B(:,i));
end

figure;plot(norms(1:2));
hold on;plot(norms(3:10));
hold on;plot(norms(11:18))
%% alignment: dot products?
%figure;plot(AnimalAlpha(:,:,2)')
%figure;plot(AnimalAlpha(:,:,5)')

nanimals = 6;
vis1_r2 = nan(nanimals,8);
vis2_r2 = nan(nanimals,8);
for animali = 1:nanimals
for i = 1:8 % with go only
    tem = corrcoef(reshape(AnimalAlpha(animali,:,1),1,[]),reshape(AnimalAlpha(animali,:,i+2),1,[]));
    vis1_r2(animali,i) = abs(tem(1,2))^1;
    tem = corrcoef(reshape(AnimalAlpha(animali,:,2),1,[]),reshape(AnimalAlpha(animali,:,i+2),1,[]));
    vis2_r2(animali,i) = abs(tem(1,2))^1;
end
end

%figure;plot(vis2_r2');%hold on;plot(vis2_r2)



vis1_r2 = nan(1,8);
vis2_r2 = nan(1,8);

for i = 1:8 % with go only
    tem = corrcoef(reshape(AnimalAlpha(:,:,1),1,[]),reshape(AnimalAlpha(:,:,i+2),1,[]));
    vis1_r2(i) = abs(tem(1,2))^1;
    tem = corrcoef(reshape(AnimalAlpha(:,:,2),1,[]),reshape(AnimalAlpha(:,:,i+2),1,[]));
    vis2_r2(i) = abs(tem(1,2))^1;
end

figure;plot(vis1_r2);hold on;plot(vis2_r2)
legend('vis1','vis2');ylabel('corcoef');xlabel('time')
%%
figure;
for i =1:10
subplot(1,10,i); shadedErrorBar([],mean(AnimalAlpha(:,:,i)),2*std(AnimalAlpha(:,:,i)/sqrt(nanimals)));ylim([0 1])
end
%% load dyn
figure;
[V,D] = eig(res.params.model.A);
real_D = abs(diag(D))';
[~,ord]=sort(real_D);
hold on;plot(real_D(ord),'k')

[V,D] = eig(res_dyn.params.model.A);
real_D = abs(diag(D))';
[~,ord]=sort(real_D);
hold on;plot(real_D(ord),'b')
legend({'constant input','changing input'})
ylabel('abs_lambda_A');xlabel('states')

figure;subplot(1,2,1);plot(res.params.model.B(:,1),'k')
subplot(1,2,2);plot(res_dyn.params.model.B(:,1:10),'b')
ylabel('B_weights');xlabel('states')

figure;subplot(1,2,1);imagesc(res_dyn.seq(45).u(1:10,:));
subplot(1,2,2);imagesc(res.seq(45).u(1,:));

%% joint

area = 'V1'; % when using fb animals (reversed stimulus is used by default)

animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66','MPV35_2'}; % remove last one?

nanimals = 6;

Q1 = [];Q2 = [];Q3 = [];Q4 = [];
for animalnum = 1:nanimals
    
    Mfile = dir(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animalnum},...
        '/Joint_trial_based_split_0/50msBins/*','.mat']);
    cd(Mfile(end).folder)
    res = load(Mfile(end-1).name);
    figure;imagesc(abs(res.params.model.A))
    title('abs(A)')
    hold on;line([0 10.5],[5.5 5.5],'color','w')
    hold on;line([5.5 5.5],[0 10.5],'color','w')
    Q1 = [Q1,abs(reshape(res.params.model.A(1:5,1:5),1,[]))];
    Q2 = [Q2,abs(reshape(res.params.model.A(1:5,6:10),1,[]))];
    Q3 = [Q3,abs(reshape(res.params.model.A(6:10,1:5),1,[]))];
    Q4 = [Q4,abs(reshape(res.params.model.A(6:10,6:10),1,[]))];
    
 %   colormap( bluewhitered(40))
end

 figure;boxplot([Q1;Q2;Q3;Q4]','notch','on','OutlierSize',0.1)