if params.model.notes.useB
  uDim = 2;
else
  uDim = 0;
end

S11 = zeros(xDim,xDim);
S01 = zeros(xDim+uDim,xDim);
S00 = zeros(xDim+uDim,xDim+uDim);

x0 = zeros(xDim,Trials);
Q0 = zeros(xDim,xDim);

Tall = [];
toSub = zeros(xDim,xDim);


uall = [];
mu0all = [];
mu1all = [];
var00 = zeros(16,16);
var01=zeros(16,16);
for tr = 1:Trials

    T = size(seq(tr).y,2);
    Tall  = [Tall T];

    if isfield(seq(tr).posterior,'Vsm')
      Vsm   = reshape(seq(tr).posterior.Vsm' ,xDim,xDim,T);
      VVsm  = reshape(seq(tr).posterior.VVsm',xDim,xDim,T-1);
    else
      Vsm   = reshape(seq(1).posterior.Vsm' ,xDim,xDim,T);
      VVsm  = reshape(seq(1).posterior.VVsm',xDim,xDim,T-1);
    end

    MUsm0 = seq(tr).posterior.xsm(:,1:T-1);
    MUsm1 = seq(tr).posterior.xsm(:,2:T);

    S11                = S11                + sum(Vsm(:,:,2:T),3)  + MUsm1*MUsm1';
    S01(1:xDim,:)      = S01(1:xDim,:)      + sum(VVsm(:,:,1:T-1),3) + MUsm0*MUsm1';
    S00(1:xDim,1:xDim) = S00(1:xDim,1:xDim) + sum(Vsm(:,:,1:T-1),3)  + MUsm0*MUsm0';

        
 
    
    if  params.model.notes.useB
      u = seq(tr).u(1:uDim,1:T-1);

            S01(1+xDim:end,:)          = S01(1+xDim:end,:)          + u*MUsm1';% masking: instead of : try 9:16 
            S00(1+xDim:end,1:xDim)     = S00(1+xDim:end,1:xDim)     + u*MUsm0';
             S00(1:xDim,1+xDim:end)     = S00(1:xDim,1+xDim:end)     + MUsm0*u';
             S00(1+xDim:end,1+xDim:end) = S00(1+xDim:end,1+xDim:end) + u*u';

      
    end
    
    
    uall = [uall,u];
    mu0all = [mu0all,MUsm0];
    mu1all = [mu1all,MUsm1];
    var00 = var00+sum(Vsm(:,:,1:T-1),3);
    var01 = var01+sum(VVsm(:,:,1:T-1),3);

    x0(:,tr) = MUsm0(:,1);
    Q0 = Q0 + Vsm(:,:,1);

end


S00 = (S00+S00')/2;
S11 = (S11+S11')/2;
%% skip
A  = S01'/S00;

B = A(:,1+xDim:end);
A = A(:,1:xDim);


B(1:8,2) = 0;

%test = (S01(1:16,1:16)'- B*u*MUsm0')/S00(1:16,1:16);
A_new = (S01(1:16,1:16)'- B*uall*mu0all')/S00(1:16,1:16);
%A = A_new
%% leave Q out for now, and 1 trial only
S01(18,1:8) = 0;
S00(18,1:8) = 0;
S00(1:8,18) = 0;


A  = S01'/S00;

B = A(:,1+xDim:end);
%B(1:8,2) = 0;

A = A(:,1:xDim);
%% derivatives

% A
difA = A*(mu0all*mu0all'+var00) - (mu0all*mu1all'+var01)' + B*uall*mu0all';

% B 
difB =A*mu0all*uall' - mu1all*uall' + B*uall*uall';

figure;subplot(1,2,1);imagesc(difA);subplot(1,2,2);imagesc(difB);
%% the optimal solution: trick: the rows are indepent in the equations above
% (although A and B are not)
A = [A1(1:8,:);A2(9:16,:)]
B = [[B1(1:8,:),zeros(8,1)];B2(9:16,:)]
% A1 from 1 input no mask, A2 from 2 inputs, second one masked
% maybe report Q_dyn here as well (?)
%%

a16=S00(1:16,1:16);
a17=S00(1:17,1:17);

a16_inv = a16^(-1);
a17_inv = a17^(-1);

%%
figure;subplot(1,2,1);
imagesc((a16 - toSub)^(-1))
subplot(1,2,2);
imagesc(a17^(-1))


%mat1 = (a16 - (MUsm0*u')*(MUsm0*u')'*((u*u')^(-1)))^(-1);
mat1 = (a16 - toSub)^(-1);
mat2 = a17^(-1);
figure;imagesc(mat1 - mat2(1:16,1:16)); % this is basically 0 (10^-18)
