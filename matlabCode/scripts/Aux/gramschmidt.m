function [U]=gramschmidt(V)

if 1
    [n,k] = size(V);
    U = zeros(n,k);
    U(:,1) = V(:,1)/norm(V(:,1));
    for i = 2:k
        U(:,i)=V(:,i);
        for j=1:i-1
            U(:,i)=U(:,i)-(U(:,j)'*U(:,i)) /(norm(U(:,j)))^2 * U(:,j);
        end
        U(:,i) = U(:,i)/norm(U(:,i));
    end
    
else
    v= V;
    k = size(v,2);
    assert(k>=2,'The input matrix must include more than one vector.');
    for ii = 1:1:k
        v(:,ii) = v(:,ii) / norm(v(:,ii));
        for jj = ii+1:1:k
            v(:,jj) = v(:,jj) - proj(v(:,ii),v(:,jj));
        end
    end
    U = v;
end


end

function w = proj(u,v)
% This function projects vector v on vector u
w = (dot(v,u) / dot(u,u)) * u;
end