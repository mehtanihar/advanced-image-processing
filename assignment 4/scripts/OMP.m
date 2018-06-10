function [ theta ] = OMP(y,A,sigma) % Plug-in m
[N,K] = (size(A));
theta= zeros(K,1);
m=size(y,1);
r=y;
T=[];
while norm(r) > sigma*sqrt(m) %while (||r||^2 > e) e ? sigma(sqrt(m)).
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

