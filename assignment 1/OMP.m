function [ theta ] = OMP(y,A) % Plug-in m
[N,K] = (size(A));
theta= zeros(K,1);
r=y;
T=[];
sigma=2;
while norm(r) > 0.000001 %while (||r||^2 > e) e ? sigma(sqrt(m)).
    g= A'*r;
    for b3=1:K
        f(b3)=abs(g(b3))/norm(A(:,b3));
    end
    
   
    jj = find(f == max(f));
    T = union(T,jj);
    theta(T)=pinv(A(:,T)) * y;  
    r=y-A*theta;  
    
end

end

