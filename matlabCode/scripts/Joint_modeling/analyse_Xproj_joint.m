animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};

%orth = 0;
figure;
%Xproj_go = cell(1,6);
%Xproj_nogo = cell(1,6);
for animali  = 1:6
    cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
    
    
    
    d = dir(['*RND*_onlyCorrect_exGrooming_go.mat']);
    rnd_vb = nan(1,length(d));
    for i = 1:length(d)
        load(d(i).name);
        vb = varBound(~isnan(varBound));
        rnd_vb(i) = vb(end);
    end
     [~,n]=max(rnd_vb);
     res = load(d(n).name);
%[~,n]=sort(rnd_vb,'descend');
%     for i = n
%         if 1%~numel(strfind(d(i).name,'onlyCorrect'))
%             res = load(d(i).name);
%             %         if max(abs(eig(res.params.model.A(1:8,1:8))))<=1
%             %             break;
%             %         end
%             [CO,W] = orth_c_svd(res.params);
%             A_new = W*res.params.model.A*(W');
%             if max(abs(eig(A_new(1:8,1:8))))<=1
%                 break
%             end
%         end
%     end
    
    % orth here
    
    A = res.params.model.A;
    
     [CO,W] = orth_c_svd(res.params);
     A_new = W*res.params.model.A*(W');
     A = A_new;
    
    [c,s,l]=pca(squeeze(A(1:8,9:16)),'Centered',0); % originally with transpose, removed the transpose dor now
    % not sure if transpose neede or not
    % hold on;plot(l/norm(l)); % almost 1d
    LMdir = c(:,1)';
    
    xproj = nan(length(res.seq),res.seq(1).T);
    for tr = 1:length(res.seq)
         if sum(sum(res.seq(tr).u(2,:),1)') == 0
            X = res.seq(tr).posterior.xsm;
            X = W * X; 
            % get only the direction of X, and only LM
%             X = X(9:16,:);
              X = X(1:8,:);
            for ti = 1:size(X,2)
                X(:,ti) = X(:,ti)./norm(X(:,ti));
            end
            
            %xproj(tr,:)=  LMdir * X(:,:);
            xproj(tr,:)= (A(1:8,9:16)* LMdir')' * X(1:8,:);
        end

    end
    subplot(1,6,animali);hold on; plot(nanmean(xproj))
    Xproj_go{animali} = xproj;
    
end
%%
figure;
for animali = 1:6
    avgo = nanmean(Xproj_go{animali});
    avnogo = nanmean(Xproj_nogo{animali});
    avgo = avgo - nanmean(avgo(20:40));
    avnogo = avnogo - nanmean(avnogo(20:40));
    subplot(1,6,animali);plot(avgo-avnogo);
end
    
    
