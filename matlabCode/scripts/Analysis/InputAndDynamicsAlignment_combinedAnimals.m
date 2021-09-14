
% alignment of lambdas (eigen values of dynamics matrix) and alphas
% (projection of input weights onto the eigen vectors of tyhe dynamics
% matrix) for all animals (in the long range area)

area = 'V1'; % when using fb animals (reversed stimulus is used by default)

animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66','MPV35_2'};

alpha = cell(1,7);
real_D = cell(1,7);

for animalnum = 1:7
    
    Mfile = dir(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animalnum},...
        '/trial_based_MS_split_1/17msBins/*_',area,'_*.mat']);
    cd(Mfile(end).folder)
    res = load(Mfile(end).name);

    [V,D] = eig(res.params.model.A);
    % do real and abs
    real_D{animalnum} = abs(diag(D))';% chnge '
    alpha{animalnum}= nan(18,size(V,2));
    for n = 1:size(V,2)
        for i = 1:18
            alpha{animalnum}(i,n) = (res.params.model.B(:,i)'/norm(res.params.model.B(:,i)))*V(:,n);
        end
    end
    alpha{animalnum} = abs(alpha{animalnum}); % real or abs
    
end

%% plots
gol = zeros(8,10);
nogol = zeros(8,10);
gov = zeros(1,10);
nogov = zeros(1,10);
figure;
for animalnum = 3%:7
    [~,ord]=sort(real_D{animalnum});
    
    subplot(2,9,9)
    plot(real_D{animalnum}(ord),'k'); 
    hold on;plot(1*(alpha{animalnum}(1,ord)),'Color',[0 1 0],'LineWidth',.01)
    ylim([-1 1])
    xlabel('state number');ylabel('visual input');
    
    subplot(2,9,18)
    plot(real_D{animalnum}(ord),'k'); 
    hold on;plot(1*(alpha{animalnum}(2,ord)),'Color',[1 0 0],'LineWidth',.01)
    ylim([-1 1])
    xlabel('state number');ylabel('visual input');
    
    gov = gov + alpha{animalnum}(1,ord);
    nogov = nogov + alpha{animalnum}(2,ord);
    for i =3:10
        subplot(2,9,i-2)
        plot(real_D{animalnum}(ord),'k');
        hold on;plot(1*(alpha{animalnum}(i,ord)),'Color',[0 1 0],'LineWidth',.01)
        ylim([-1 1])
        gol(i-2,:) = gol(i-2,:) + alpha{animalnum}(i,ord);
        xlabel('state number');%ylabel('feedback input');       
    end
    for i =11:18
        subplot(2,9,i-1)
        plot(real_D{animalnum}(ord),'k');
        hold on;plot(1*(alpha{animalnum}(i,ord)),'Color',[1 0 0],'LineWidth',.01)
        ylim([-1 1])
        nogol(i-10,:) = nogol(i-10,:) + alpha{animalnum}(i,ord);
        xlabel('state number');%ylabel('feedback input');
    end
end
gol = gol/7;
nogol = nogol/7;
gov = gov/7;
nogov = nogov/7;

for i = 1:8
    subplot(2,9,i); plot(gol(i,:),'Color',[0 0.8 0],'LineWidth',3)
    subplot(2,9,i+9); plot(nogol(i,:),'Color',[0.8 0 0],'LineWidth',3)
end

%%% visual stimulus now:

subplot(2,9,9); plot(gov,'Color',[0 0.8 0],'LineWidth',3)
subplot(2,9,18); plot(nogov,'Color',[0.8 0 0],'LineWidth',3)


