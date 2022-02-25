
% alignment of lambdas (eigen values of dynamics matrix) and alphas
% (projection of input weights onto the eigen vectors of tyhe dynamics
% matrix) for all animals (in the long range area)

area = 'V1'; % when using fb animals (reversed stimulus is used by default)

animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66','MPV35_2'}; % remove last one?

nanimals = 1;

alpha = cell(1,nanimals);
real_D = cell(1,nanimals);

AllNorms = nan(nanimals,18); % only works for when no dynInput

for animalnum = 1:nanimals
    
%     Mfile = dir(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animalnum},...
%         '/trial_based_MS_split_1/17msBins/*_',area,'_nState8*.mat']);
         Mfile = dir(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animalnum},...
             '/trial_based_split_1/17msBins/*','_DiscInput_onlyCorrect_exGrooming_RND4.mat']); % 17ms BaselineU = 0
%         Mfile = dir(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animalnum},...
%             '/trial_based_split_1/50msBins/*','02_07_*.mat']);
    cd(Mfile(end).folder)
    res = load(Mfile(end).name);

    [V,D] = eig(res.params.model.A);
    % do real and abs
    real_D{animalnum} = abs(diag(D))';% chnge '
%     alpha{animalnum}= nan(18,size(V,2));
%     for n = 1:size(V,2)
%         for i = 1:18
%             alpha{animalnum}(i,n) = (res.params.model.B(:,i)'/norm(res.params.model.B(:,i)))*V(:,n);
%         end
%     end
        B = res.params.model.B;%nan(size(res.params.model.B));
        size(B,2)
%         for i = 1:size(B,2)
%               B(:,i) = (res.params.model.B(:,i)/norm(res.params.model.B(:,i)));
%         end
       
%      alpha{animalnum}= nan(18,size(V,2));
%      for n = 1:size(V,2)
%          for i = 1:18
%              alpha{animalnum}(i,n) = (res.params.model.B(:,i)'/norm(res.params.model.B(:,i)))*V(:,n);
%          end
%      end

    alpha{animalnum}= (inv(V)*B)';%
    alpha{animalnum} = abs(alpha{animalnum}); % real or abs
    
     % forDynInput - 50 ms -- res end-1
 %   alpha{animalnum} = [alpha{animalnum}(1,:);alpha{animalnum}(11,:);alpha{animalnum}(end-15:end,:)];
    
    for i =1:size(alpha{animalnum},1)
        AllNorms(animalnum,i) = norm(alpha{animalnum}(i,:));
        alpha{animalnum}(i,:) = alpha{animalnum}(i,:)/norm(alpha{animalnum}(i,:));
        
    end
    
   
end

% Todo: normalize per animals
yl = 1; % 0.05 or 2
figure;plot(AllNorms');
%figure;shadedErrorBar([],mean(AllNorms(:,3:10)),2*std(AllNorms(:,3:10))/sqrt(nanimals),'k',0)
%xlabel('time');ylabel('norm_B_go, mean ci')
%%% intermediate plots
figure
for animalnum = 1:nanimals
    [~,ord]=sort(real_D{animalnum});
 hold on;plot(real_D{animalnum}(ord),'k')
 xlabel('state number')
 ylabel('abs(lambda)')
end
%%% plots
nstates = 8;

gol = zeros(8,nstates);
nogol = zeros(8,nstates);
gov = zeros(1,nstates);
nogov = zeros(1,nstates);
figure;

AnimalAlpha = nan(nanimals,nstates,18);

for animalnum = 1:nanimals
    [~,ord]=sort(real_D{animalnum});
    
    s = subplot(2,9,9);
    plot(real_D{animalnum}(ord),'k'); 
    hold on;plot(1*(alpha{animalnum}(1,ord)),'Color',[0.5 0.5 0.5],'LineWidth',.01)
    AnimalAlpha(animalnum,:,:) = alpha{animalnum}(:,ord)';
    ylim([0 yl])
    xlabel('state number');
    s.YTick = [0 1];
    s.YTickLabel = [0 1];
    s.YTick = [1 8];
    s.YTickLabel = [ 1 8];   
    ylabel('abs(alpha)');    
    s.Title.String = ('visual input');
    s.Title.FontWeight = 'normal';
    s.Title.Position(2) = s.Title.Position(2)+0.005;
    
    s=subplot(2,9,18);
    plot(real_D{animalnum}(ord),'k'); 
    hold on;plot(1*(alpha{animalnum}(2,ord)),'Color',[.5 0.5 0.5],'LineWidth',.01)
    ylim([0 yl])
    xlabel('state number');
    s.YTick = [0 1];
    s.YTickLabel = [0 1];
    s.YTick = [1 8];
    s.YTickLabel = [ 1 8];   
    ylabel('abs(alpha)');    
    s.Title.String = ('visual input');
    s.Title.FontWeight = 'normal';
    s.Title.Position(2) = s.Title.Position(2)+0.005;
    
    gov = gov + alpha{animalnum}(1,ord);
    nogov = nogov + alpha{animalnum}(2,ord);
    for i =3:10
        s = subplot(2,9,i-2);
        plot(real_D{animalnum}(ord),'k');
        hold on;plot(1*(alpha{animalnum}(i,ord)),'Color',[0.5 .5 0.5],'LineWidth',.01)
        ylim([0 yl])
        gol(i-2,:) = gol(i-2,:) + alpha{animalnum}(i,ord);
        xlabel('state number');%ylabel('feedback input');   
        
        s.YTick = [0 1];
        s.YTickLabel = [0 1];
        s.YTick = [1 8];
        s.YTickLabel = [ 1 8];
        if i ==3
        ylabel('abs(alpha)');
        end
        s.Title.String = (['delay',num2str(i-2)]);
        s.Title.FontWeight = 'normal';
        s.Title.Position(2) = s.Title.Position(2)+0.01;
    end
    for i =11:18
        s= subplot(2,9,i-1);
        plot(real_D{animalnum}(ord),'k');
        hold on;plot(1*(alpha{animalnum}(i,ord)),'Color',[.5 0.5 0.5],'LineWidth',.01)
        ylim([0 yl])
        nogol(i-10,:) = nogol(i-10,:) + alpha{animalnum}(i,ord);
        xlabel('state number');%ylabel('feedback input');
        
        s.YTick = [0 1];
        s.YTickLabel = [0 1];
        s.YTick = [1 8];
        s.YTickLabel = [ 1 8];
        if i == 11
        ylabel('abs(alpha)');
        end
        s.Title.String = (['delay',num2str(i-10)]);
        s.Title.FontWeight = 'normal';
        s.Title.Position(2) = s.Title.Position(2)+0.01;

     end
 end
gol = gol/nanimals;
nogol = nogol/nanimals;
gov = gov/nanimals;
nogov = nogov/nanimals;

% Add anovas
for i = 1:8
subplot(2,9,i); hold on; text(0,.9*yl,num2str(anova1(AnimalAlpha(:,:,2+i),[],'off')));
subplot(2,9,i+9); hold on; text(0,.9*yl,num2str(anova1(AnimalAlpha(:,:,10+i),[],'off')));
end
subplot(2,9,9); hold on; text(0,.9*yl,num2str(anova1(AnimalAlpha(:,:,1),[],'off')));
subplot(2,9,18); hold on; text(0,.9*yl,num2str(anova1(AnimalAlpha(:,:,2),[],'off')));


for i = 1:8
    subplot(2,9,i); plot(gol(i,:),'Color',[0 0.8 0],'LineWidth',3)
    subplot(2,9,i+9); plot(nogol(i,:),'Color',[0.8 0 0],'LineWidth',3)
end

%%% visual stimulus now:

subplot(2,9,9); plot(gov,'Color',[0 0 0.8],'LineWidth',3)
subplot(2,9,18); plot(nogov,'Color',[0 0 0.8],'LineWidth',3)

%%
% figure;scatter(cell2mat(cellfun(@(x) x(2,ord),alpha,'UniformOutput',0)),...
%     cell2mat(cellfun(@(x) x(13,ord),alpha,'UniformOutput',0)))
% probably not corect cuz the discr input is ifferent directions for go vs
% nogo
cog = [];
cong = [];
for i = 1:8
    ccfg = corrcoef(cell2mat(cellfun(@(x) x(1,:),alpha,'UniformOutput',0)),...
    cell2mat(cellfun(@(x) x(2+i,:),alpha,'UniformOutput',0)));%2 or 10

    ccfng = corrcoef(cell2mat(cellfun(@(x) x(1,:),alpha,'UniformOutput',0)),...
    cell2mat(cellfun(@(x) x(10+i,:),alpha,'UniformOutput',0)));%2 or 10
    
    cog(i) = (ccfg(1,2)+ ccfng(1,2))/2;
    
    ccfg = corrcoef(cell2mat(cellfun(@(x) x(2,:),alpha,'UniformOutput',0)),...
    cell2mat(cellfun(@(x) x(2+i,:),alpha,'UniformOutput',0)));
    ccfng = corrcoef(cell2mat(cellfun(@(x) x(2,:),alpha,'UniformOutput',0)),...
    cell2mat(cellfun(@(x) x(10+i,:),alpha,'UniformOutput',0)));
    cong(i) = (ccfg(1,2)+ ccfng(1,2))/2;
end

figure;plot(cog,'g');hold on;plot(cong,'r')
% add corcoef

%%
allLambda = cell2mat(real_D);
allAlpha = cell2mat(cellfun(@(x) x(5,:),alpha,'UniformOutput',0));
corrcoef(allLambda,allAlpha)
figure;scatter(allLambda,allAlpha)
%% here the best is Disc_dyn_ and the first 3 animals. so maybe exgrooming? but only so if vis 2 is sum not dif. which one makes sense? 
nanimals =1;
vis1 = nan(nanimals,8);
vis2 = nan(nanimals,8);
for animalnum = 1:nanimals
    
        Mfile = dir(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animalnum},...
            '/trial_based_split_1/17msBins/*','DiscInput_exGrooming.mat']);
    cd(Mfile(end).folder)
    res = load(Mfile(end).name);
    
%             Mfile = dir(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animalnum},...
%             '/trial_based_split_1/50msBins/*','02_07_*.mat']);
%     cd(Mfile(end).folder)
%     res = load(Mfile(end-1).name);

%     [V,D] = eig(res.params.model.A);
%     % do real and abs
%     real_D{animalnum} = abs(diag(D))';% chnge '
%     alpha{animalnum}= nan(18,size(V,2));
%     for n = 1:size(V,2)
%         for i = 1:18
%             alpha{animalnum}(i,n) = (res.params.model.B(:,i)'/norm(res.params.model.B(:,i)))*V(:,n);
%         end
%     end
        B = res.params.model.B;%nan(size(res.params.model.B));
        for i = 1:size(B,2)
              B(:,i) = (res.params.model.B(:,i)/norm(res.params.model.B(:,i)));
        end
        
        % go trials only
%         vis1(animalnum,:) = (B(:,1)'*B(:,3:10));
%         vis2(animalnum,:) = (B(:,2)'*B(:,3:10));
%         
%          vis1(animalnum,:) = (B(:,1)'*B(:,11:18));
%         vis2(animalnum,:) = (-B(:,2)'*B(:,11:18));


%         vis1(animalnum,:) = (B(:,1)'*B(:,21:28));
%         vis2(animalnum,:) = (B(:,11)'*B(:,21:28));
        
%           vis1(animalnum,:) = (B(:,1)'*B(:,29:36));
%          vis2(animalnum,:) = (B(:,11)'*B(:,29:36)); %- needed?
         
%          
%           vis1(animalnum,:) = (B(:,1)'*B(:,29:36)) +  (B(:,1)'*B(:,21:28));
%           vis2(animalnum,:) = (B(:,11)'*B(:,29:36)) + (B(:,11)'*B(:,21:28));

             vis1(animalnum,:) = (B(:,1)'*B(:,3:10));% +  (B(:,1)'*B(:,11:18));
             vis2(animalnum,:) = (B(:,2)'*B(:,3:10));% - (B(:,2)'*B(:,11:18));
        
end
% for i = 1:size(vis1,1)
%     vis1(i,:) = vis1(i,:)/norm(vis1(i,:));
%     vis2(i,:) = vis2(i,:)/norm(vis2(i,:));
% end
figure;subplot(1,2,1);plot(vis1','k');hold on;plot(mean(vis1),'k','LineWidth',2)
subplot(1,2,2);plot(vis2','k');xlabel('time');hold on;plot(mean(vis2),'k','LineWidth',2)
%%
a = nan(nanimals,8); % 8 is ntimewindows
for animalnum = 1:nanimals
    [~,ord]=sort(real_D{animalnum});
    a(animalnum,:)=alpha{animalnum}(3:10,ord(5))';
end
figure;plot(a')
