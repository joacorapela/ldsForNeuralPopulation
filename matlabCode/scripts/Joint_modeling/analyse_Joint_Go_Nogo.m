%%% At the moment, for each animal, rand0 is selected. To be added later:
%%% going through all rounds and picking the best model

%% visualizing simgle cell responses - the only part requiring seq
% load one of the go or nogo files
% first panel is cntrl trials, second is laser at lag 2 - [the selectoin of
% laser only works for 17 ms bins].
figure;
for i = 1:10
    heldoutN =randi(size(params.model.C,1))
    pred_bs = nan(length(seq),seq(1).T); % trials*timebins
    yOrig_bs = nan(length(seq),seq(1).T); % trials*timebins
    pred_ls = nan(length(seq),seq(1).T); % trials*timebins
    yOrig_ls = nan(length(seq),seq(1).T); % trials*timebins
    
    uv_bs = nan(length(seq),seq(1).T);
    ul_bs = nan(length(seq),seq(1).T);
    uv_ls = nan(length(seq),seq(1).T);
    ul_ls = nan(length(seq),seq(1).T);
    for tr = 1:length(seq)
        % for trial tr:
        z = params.model.C(heldoutN,:) *seq(tr).posterior.xsm + ...
            params.model.d(heldoutN);
        
        %  if sum(sum(seq(tr).u(2,:),1)) == 0
        pred_bs(tr,:) = exp(z);
        yOrig_bs(tr,:) = seq(tr).y(heldoutN,:);
        uv_bs(tr,:) = seq(tr).u(1,:);
        %        ul_bs(tr,:) = seq(tr).u(2,:);
        %     elseif(seq(tr).u(2,62) - seq(tr).u(2,61)) > 0 % change to else, if you want all silencing onsets not only lag 2
        %         pred_ls(tr,:) = exp(z);
        %         yOrig_ls(tr,:) = seq(tr).y(heldoutN,:);
        %         uv_ls(tr,:) = seq(tr).u(1,:);
        %         ul_ls(tr,:) = seq(tr).u(2,:);
        %     end
    end
    
    
    subplot(10,1,i);plot(nanmean(pred_bs,1));hold on;plot(nanmean(yOrig_bs,1));
    hold on;plot(nanmean(uv_bs,1));hold on;plot(nanmean(ul_bs,1));legend({'model prediction','original trace', 'visual input'});
end
figure;
subplot(1,2,1);plot(nanmean(pred_bs,1));hold on;plot(nanmean(yOrig_bs,1));
hold on;plot(nanmean(uv_bs,1));hold on;plot(nanmean(ul_bs,1));legend({'bs prediction','bs original trace', 'vis', 'laser'});
subplot(1,2,2);plot(nanmean(pred_ls,1));hold on;plot(nanmean(yOrig_ls,1));
hold on;plot(nanmean(uv_ls,1));hold on;plot(nanmean(ul_ls,1));legend({'ls prediction','ls original trace', 'vis', 'laser'})
%% fit examples
cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/','MPV17','/Joint_trial_based_splitContext_CntrlOnly/17msBins/'])
res = load('Joint_PLDSfitRes_22_02_20_06_45_27_RND4_onlyCorrect_exGrooming_go.mat');
orth = 2;

if orth == 2
    [CO,W] = orth_c_svd(res.params);
end
figure;

rng(3,'twister')
for i = 1:10
    
    heldoutN =randi(size(res.params.model.C,1))
    
    pred_bs = nan(length(res.seq),res.seq(1).T); % trials*timebins
    yOrig_bs = nan(length(res.seq),res.seq(1).T); % trials*timebins
    
    uv_bs = nan(length(res.seq),res.seq(1).T);
    for tr = 1:length(res.seq)
        % for trial tr:
        z = res.params.model.C(heldoutN,:) *res.seq(tr).posterior.xsm + ...
            res.params.model.d(heldoutN);
        
        %         z = CO(heldoutN,:)*W*res.seq(tr).posterior.xsm + ...
        %         res.params.model.d(heldoutN);
        
        pred_bs(tr,:) = exp(z);
        yOrig_bs(tr,:) = res.seq(tr).y(heldoutN,:);
        uv_bs(tr,:) = res.seq(tr).u(1,:);
        
    end
    
    
    subplot(10,1,i);plot(nanmean(pred_bs,1));hold on;plot(nanmean(yOrig_bs,1));
    hold on;plot(nanmean(uv_bs,1));hold on;plot(nanmean(ul_bs,1));legend({'model prediction','original trace', 'visual input'},'Location','northwest');
end
%% loading params
animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};

orth = 2;

allA = nan(6,16,16,2);
allB = nan(6,16,2);
for animali  = 1:6
    cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
    disp('-----')
    
    
     d = dir(['*RND*_onlyCorrect_exGrooming_go.mat']);
  %  d = dir(['*RND*_exGrooming_go.mat']);
    
    rnd_vb = nan(1,length(d));
    for i = 1:length(d)
        if  numel(strfind(d(i).name,'onlyCorrect'))
        load(d(i).name);
        vb = varBound(~isnan(varBound));
        rnd_vb(i) = vb(end);
        end
    end
%              [~,n]=max(rnd_vb);
%              res = load(d(n).name);
    [~,n]=sort(rnd_vb,'descend');
    for i = n
        if 1%~numel(strfind(d(i).name,'onlyCorrect'))
            i
            res = load(d(i).name);
            %                 if max(abs(eig(res.params.model.A(1:8,1:8))))<=1
            %                     break;
            %                 end
            
            [CO,W] = orth_c_svd(res.params);
            A_new = W*res.params.model.A*(W');
            if max(abs(eig(A_new(1:8,1:8))))<=1;%+0.001
                break
            end
        end
    end
    if orth ==1
        [CO]=gramschmidt(res.params.model.C); % check this later -- importantly CO still keeps V1 and LM separate - maybe check this for all animals
        W = CO'*res.params.model.C;
        A_new = W*res.params.model.A*inv(W);
        allA(animali,:,:,1) = A_new;%params.model.A;
    elseif orth == 2
        [CO,W] = orth_c_svd(res.params);
        A_new = W*res.params.model.A*(W');
        allA(animali,:,:,1) = A_new;
        allB(animali,:,1) =  W*res.params.model.B(:,1);
    else
        allA(animali,:,:,1) = res.params.model.A;
        allB(animali,:,1) = res.params.model.B(:,1);
    end
    
    
    d = dir(['*RND*_onlyCorrect_exGrooming_nogo.mat']);
   % d = dir(['*RND*_exGrooming_nogo.mat']);
    
    rnd_vb = nan(1,length(d));
    for i = 1:length(d)
         if  numel(strfind(d(i).name,'onlyCorrect'))
        load(d(i).name);
        vb = varBound(~isnan(varBound));
         end
        rnd_vb(i) = vb(end);
    end
%              [~,n]=max(rnd_vb);
%              res=load(d(n).name);
    [~,n]=sort(rnd_vb,'descend');
    for i = n
        if 1%~numel(strfind(d(i).name,'onlyCorrect'))
            res = load(d(i).name);
            %         if max(abs(eig(res.params.model.A(1:8,1:8))))<=1
            %             break;
            %         end
            [CO,W] = orth_c_svd(res.params);
            A_new = W*res.params.model.A*(W');
            if max(abs(eig(A_new(1:8,1:8))))<=1
                break
            end
        end
    end
    
    if orth == 1
        [CO]=gramschmidt(res.params.model.C); % check this later -- importantly CO still keeps V1 and LM separate - maybe check this for all animals
        W = CO'*res.params.model.C;
        A_new = W*res.params.model.A*inv(W);
        allA(animali,:,:,2) = A_new;%params.model.A;
    elseif orth == 2
        [CO,W] = orth_c_svd(res.params);
        A_new = W*res.params.model.A*(W');
        allA(animali,:,:,2) = A_new;
        allB(animali,:,2) =  W*res.params.model.B(:,1);
    else
        allA(animali,:,:,2) = res.params.model.A;
        allB(animali,:,2) =  res.params.model.B(:,1);
    end
    
end


%%% to test orth:
% figure;imagesc(params.model.C - CO*W)
%% visuallization of A
animali = 6
figure;subplot(1,2,1);imagesc(squeeze(allA(animali,:,:,1)));caxis([-0.5 0.5])
subplot(1,2,2);imagesc(squeeze(allA(animali,:,:,2)));caxis([-0.5 0.5])


%%

wgo=[reshape(abs(allA(:,1:8,1:8,1)),1,[]);reshape(abs(allA(:,1:8,9:16,1)),1,[]);reshape(abs(allA(:,9:16,1:8,1)),1,[]);reshape(abs(allA(:,9:16,9:16,1)),1,[])];
wnogo=[reshape(abs(allA(:,1:8,1:8,2)),1,[]);reshape(abs(allA(:,1:8,9:16,2)),1,[]);reshape(abs(allA(:,9:16,1:8,2)),1,[]);reshape(abs(allA(:,9:16,9:16,2)),1,[])];

figure;subplot(1,2,1);boxplot(wgo');subplot(1,2,2);boxplot(wnogo')
ranksum(wgo(2,:),wnogo(2,:))
%% dimensionality of interareal and intraareal subspaces: both go and nogo pooled



varExpPC_LV = nan(2,6,8);
varExpPC_LL = nan(2,6,8);
varExpPC_VV = nan(2,6,8);
for animali  =[1 2 3 4 5 6]%1:6
    for lag = 1:2      
        [c,s,l]=pca(squeeze(allA(animali,1:8,9:16,lag)),'Centered',0); % originally it was transpose, makes sense?
        varExpPC_LV(lag,animali,:) = l/sum(l); % almost 1d        
        
        [c,s,l]=pca(squeeze(allA(animali,1:8,1:8,lag)),'Centered',0); % originally it was transpose, makes sense?
        varExpPC_VV(lag,animali,:) = l/sum(l); % almost 1d        
        
        [c,s,l]=pca(squeeze(allA(animali,9:16,9:16,lag)),'Centered',0); % originally it was transpose, makes sense?
        varExpPC_LL(lag,animali,:) = l/sum(l); % almost 1d        
    end
end

f=figure;boxplot([reshape(varExpPC_LV(:,:,1),1,[]);reshape(varExpPC_LV(:,:,2),1,[]);reshape(varExpPC_LV(:,:,3),1,[]);reshape(varExpPC_LV(:,:,4),1,[]);reshape(varExpPC_LV(:,:,5),1,[]);...
    reshape(varExpPC_LV(:,:,6),1,[]);reshape(varExpPC_LV(:,:,7),1,[]);reshape(varExpPC_LV(:,:,8),1,[])]','Widths',.4);
set(f.Children,'YLim',[0,1]);f.Children.Title.String = 'A_{lv}';
f.Children.YLabel.String = 'variance explained(lambda_i)'; f.Children.XLabel.String = 'i';set(f,'color','white');

f=figure;boxplot([reshape(varExpPC_VV(:,:,1),1,[]);reshape(varExpPC_VV(:,:,2),1,[]);reshape(varExpPC_VV(:,:,3),1,[]);reshape(varExpPC_VV(:,:,4),1,[]);reshape(varExpPC_VV(:,:,5),1,[]);...
    reshape(varExpPC_VV(:,:,6),1,[]);reshape(varExpPC_VV(:,:,7),1,[]);reshape(varExpPC_VV(:,:,8),1,[])]','Widths',.4);
set(f.Children,'YLim',[0,.5]);f.Children.Title.String = 'A_{vv}';
f.Children.YLabel.String = 'variance explained(lambda_i)'; f.Children.XLabel.String = 'i';set(f,'color','white');

f=figure;boxplot([reshape(varExpPC_LL(:,:,1),1,[]);reshape(varExpPC_LL(:,:,2),1,[]);reshape(varExpPC_LL(:,:,3),1,[]);reshape(varExpPC_LL(:,:,4),1,[]);reshape(varExpPC_LL(:,:,5),1,[]);...
    reshape(varExpPC_LL(:,:,6),1,[]);reshape(varExpPC_LL(:,:,7),1,[]);reshape(varExpPC_LL(:,:,8),1,[])]','Widths',.4);
set(f.Children,'YLim',[0,.5]);f.Children.Title.String = 'A_{ll}';
f.Children.YLabel.String = 'variance explained(lambda_i)'; f.Children.XLabel.String = 'i';set(f,'color','white');

%% determine go or nogo:

lag =1 ; % 1=go, 2=nogo
type = 'box2';%'box or 'erbar'

mat = squeeze(varExpPC_LV(lag,:,:));
f=figure;
if strcmp(type,'box')
    boxplot(mat,'Widths',.4);
else
   errorbar([],nanmean(mat,1),2*nanstd(mat)/sqrt(size(mat,1)),'color','k')
end
set(f.Children,'YLim',[0,1]);set(f.Children,'XLim',[0,9]);f.Children.Title.String = 'A_{lv}';
f.Children.YLabel.String = 'variance explained(lambda_i)'; f.Children.XLabel.String = 'i';set(f,'color','white');

mat = squeeze(varExpPC_VV(lag,:,:));
f=figure;
if strcmp(type,'box')
    boxplot(mat,'Widths',.4);
else
   errorbar([],nanmean(mat,1),2*nanstd(mat)/sqrt(size(mat,1)),'color','k')
end
set(f.Children,'YLim',[0,.5]);set(f.Children,'XLim',[0,9]);f.Children.Title.String = 'A_{vv}';
f.Children.YLabel.String = 'variance explained(lambda_i)'; f.Children.XLabel.String = 'i';set(f,'color','white');


mat = squeeze(varExpPC_LL(lag,:,:));
f=figure;
if strcmp(type,'box')
    boxplot(mat,'Widths',.4);
else
   errorbar([],nanmean(mat,1),2*nanstd(mat)/sqrt(size(mat,1)),'color','k')
end
set(f.Children,'YLim',[0,.5]);set(f.Children,'XLim',[0,9]);f.Children.Title.String = 'A_{ll}';
f.Children.YLabel.String = 'variance explained(lambda_i)'; f.Children.XLabel.String = 'i';set(f,'color','white');


%% go vs nogo comparison
f = figure; set(f,'color','white');
ComMat = varExpPC_VV;

for lag = 1:2
    mat = squeeze(ComMat(lag,:,:));
    if lag == 1
   % hold on;plot(nanmean(mat,1),'g');
    hold on;errorbar((1:8)+0.1,nanmean(mat,1),2*nanstd(mat)/sqrt(size(mat,1)),'color','g')
   % hold on;boxplot(mat,'Widths',.4,'Colors','g','Positions',(1:8)-0.2);

    else
       %  hold on;boxplot(mat,'Widths',.4,'Colors','r','Positions',(1:8)+0.2);
         hold on;errorbar((1:8)-0.1,nanmean(mat,1),2*nanstd(mat)/sqrt(size(mat,1)),'color','r')
    end
    
end
for i = 1%:8
    p = signrank(ComMat(1,:,i),ComMat(2,:,i));
    hold on;text(i,0.75,['p = ',num2str(p)])
end

%% V1 local modes and interareal
%here we see the mode that is affected in V1 is aligned with fat or slow
%local direction. Important part is go/nogo difference: if there is a
%consistent pattern in go and not nogo

f = figure;
all_proj = [];
all_lambda_mag = [];
all_lambda_ang = [];
vis_vs_sharedMode = nan(2,6);
av = nan(1,8);

norm_B = nan(2,6); % have to show that the lack of pattern in nogo is not due to noise and B being small
for animali  =[1 2 3 4 5 6]%1:6
    
    for lag = 1:2
        % LM pca
        
        %         [c,s,l]=pca(squeeze(allA(animali,1:8,1:8,lag))');
        %         LMdir = c(:,1);
        %         B = squeeze(allA(animali,1:8,1:8,lag))*LMdir; % I think this is the direction in V1 that things mve along, when LM activity moves along the com dimension
        
        
        [c,s,l]=pca(squeeze(allA(animali,1:8,9:16,lag)),'Centered',0); % originally it was transpose, makes sense?
        % not sure if transpose neede or not
        % subplot(2,6,animali);hold on;plot(l/norm(l)); % almost 1d
        LMdir = c(:,1);
        
        B = s(:,1);% or squeeze(allA(animali,1:8,9:16,lag))*LMdir; % I think this is the direction in V1 that things mve along, when LM activity moves along the com dimension
        vis_vs_sharedMode(lag,animali) = (B/norm(B))'*(allB(animali,1:8,lag)'/norm(allB(animali,1:8,lag)'));
        
        norm_B(lag,animali) = norm(B);
        B = B/norm(B);
        
        [u,v] = eig(squeeze(allA(animali,1:8,1:8,lag)));
        [~,Ord]=sort(abs((diag(v))));
        proj = abs(inv(u(:,:)) * B);
        proj = proj(Ord);
      % subplot(2,6,animali)
        if lag == 1
                subplot(1,2,1)
            hold on;plot(proj./norm(proj),'color',[.7 .7 .7]) % or green
            hold on;plot(abs((diag(v(Ord,Ord)))),'k')
        else
                 subplot(1,2,2)
            hold on;plot(proj./norm(proj),'r')
            hold on;plot(abs((diag(v(Ord,Ord)))),'k')
        end
        
        
        % subplot(2,6,animali+6); hold on;scatter(abs(diag(v(Ord,Ord))),proj./norm(proj)); hold on;plot(abs(diag(v(Ord,Ord))),proj./norm(proj))
        
        if lag ==1
            av = [av;(proj./norm(proj))'];
            all_proj = [all_proj;proj./norm(proj)];
            all_lambda_mag = [all_lambda_mag;abs((diag(v(Ord,Ord))))];
            all_lambda_ang = [all_lambda_ang;imag((diag(v(Ord,Ord))))];
            
        end
        
    end
    legend
end

av(1,:) = [];
s = subplot(1,2,1); %hold on;plot(nanmean(av),'Color',[0 .7 0],'LineWidth',2);
s.YLabel.String = '|lambda|'; s.XLabel.String = 'states';set(f,'color','white');

f = figure; scatter(all_lambda_mag,all_proj,100,'.')
f.Children.YLabel.String = 'alpha'; f.Children.XLabel.String = '|lambda|';set(f,'color','white');
[r,p]= corrcoef(all_lambda_mag,all_proj)
%% LM local modes and interareal

figure;
all_proj = [];
all_lambda_mag = [];
all_lambda_ang = [];
vis_vs_sharedMode = nan(2,6);

varExpPC_LV = nan(2,6,8);


for animali  = [1 2 3 4 5 6]%1:6
    subplot(1,6,animali)
    for lag = 1:2
        % LM pca
        
        [c,s,l]=pca(squeeze(allA(animali,1:8,9:16,lag)),'Centered',0);
        % not sure if transpose neede or not
      %   hold on;plot(l/norm(l)); % almost 1d
        varExpPC_LV(lag,animali,:) = l/sum(l);
        LMdir = c(:,1);
        vis_vs_sharedMode(lag,animali) = (LMdir/norm(LMdir))'*(allB(animali,9:16,lag)'/norm(allB(animali,9:16,lag)'));
        
        B = LMdir;
        %B = LMdir;
         B = B/norm(B);
         
        [u,v] = eig(squeeze(allA(animali,9:16,9:16,lag)));
        [~,Ord]=sort(abs((diag(v))));
        proj = abs(inv(u(:,Ord)) * B);
         hold on;plot(proj./norm(proj))
         hold on;plot(abs((diag(v(Ord,Ord)))),'k')
         
        if lag ==1
            all_proj = [all_proj;proj./norm(proj)];
            all_lambda_mag = [all_lambda_mag;abs((diag(v(Ord,Ord))))];
            all_lambda_ang = [all_lambda_ang;imag((diag(v(Ord,Ord))))];
            
        end
        
    end
    legend
end


figure;scatter(all_lambda_mag,all_proj)
[r,p]= corrcoef(all_lambda_mag,all_proj)

figure;boxplot([reshape(varExpPC_LV(:,:,1),1,[]);reshape(varExpPC_LV(:,:,2),1,[]);reshape(varExpPC_LV(:,:,3),1,[]);reshape(varExpPC_LV(:,:,4),1,[]);reshape(varExpPC_LV(:,:,5),1,[]);...
    reshape(varExpPC_LV(:,:,6),1,[]);reshape(varExpPC_LV(:,:,7),1,[]);reshape(varExpPC_LV(:,:,8),1,[])]')

%% maybe check alignment with inputs as well