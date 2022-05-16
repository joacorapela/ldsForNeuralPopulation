% calculating Q_dyn (terms depending on A,B,Q)

% marginal likelihood p({x},{y}), the terms that involve A,B,Q

res = out_c_m2;
T = res.seq(1).T;
ML = 0;

for k = 1:numel(res.seq)
    
    x0 = res.seq(k).posterior.xsm(:,1:T-1);
    x1 = res.seq(k).posterior.xsm(:,2:T);
    u0 = res.seq(k).u(:,1:T-1);
    ML = ML -.5*trace((x1 - (res.params.model.A*x0+res.params.model.B*u0))'*((res.params.model.Q)^(-1))*(x1 - (res.params.model.A*x0+res.params.model.B*u0))) - ((T-1)/2)*log(det(res.params.model.Q));
end
ML

% trace is equivalnt to this:
% s = 0;
% for i =1:T-1
%     s=s+((x1(:,i) - res.params.model.A*x0(:,i))'*((res.params.model.Q)^(-1))*(x1(:,i) - res.params.model.A*x0(:,i)));
% end
% -0.5*s