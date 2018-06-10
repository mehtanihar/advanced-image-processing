%% PART A

K=20;
p=100;
%Generate D from Normal distribution
D=randn(p,K);
%Unit normalize columns of D
for i=1:K
   D(:,i)=D(:,i)/norm(D(:,i));
end
N=100;
x=zeros(p,N);
for i=1:N
    %Choose random 5 columns of D
    col_indices=randi([1,K],5,1);
    for j=1:5
       x(:,i)=x(:,i)+rand()*10*D(:,col_indices(j)); 
    end
end

m=10; %range of values of m=[10,20,30,40,50,70,90]
mu=0;
f=0.001; %range of values of m=[0.001,0.01,0.02,0.05,0.1,0.3]
sigma=0;
phi=(randi([0,1],[N,m,p])-0.5)/(0.5*sqrt(m));

for i=1:N
    sigma=sigma+sum(abs(reshape(phi(i,:,:),[m,p])*x(:,i)));
end
sigma=f*sigma/(m*N);
y=zeros(m,N);
noise=normrnd(mu,sigma,m,N);
for i=1:N
    y(:,i)=reshape(phi(i,:,:),[m,p])*x(:,i)+noise(:,i);
end

A=zeros(N,m,K);
D_prev=zeros(p,K);
iter=0;
x_prev=zeros(p,N);
x_pred=x;

while(norm(x_prev-x_pred)>0.01)
    iter=iter+1;
    disp(norm(x_prev-x_pred));
    x_prev=x_pred;
    %Step 1: Sparse Coding
    for i=1:N
       A(i,:,:)=reshape(phi(i,:,:),[m,p])*D;
    end
    S=zeros(K,N);
    for i=1:N
        S(:,i)=OMP(y(:,i),reshape(A(i,:,:),[m,K]),0.00001);
    end
    

    %Step 2: Dictionary update
    for k=1:K
        
        G=zeros(p,p);
        h=zeros(p,1);
        indices=find(S(k,:)~=0);
       
        for it=1:length(indices)
           i=indices(it);
           phi_i=reshape(phi(i,:,:),[m,p]);
           y_k=0;
           for k1=1:K
               if(k1~=k)
                    y_k=y_k+D(:,k1)*S(k1,i);
               end
           end
           y_k=y(:,i)-phi_i*y_k;
           G = G + phi_i'*phi_i*(S(k,i)^2);
           h=h+(phi_i'*y_k*S(k,i));
        end

        dk=pinv(G)*h;
        
        dk=dk/norm(dk);
       
        for it=1:length(indices)
           i=indices(it);
           phi_i=reshape(phi(i,:,:),[m,p]);
           S(k,i)=(dk'*phi_i'*y_k)/(dk'*(phi_i'*phi_i)*dk);
                  
        end
        D(:,k)=dk;
        
    end
    x_pred=D*S;
    
    
end
%% 

% 
% 
% A=D;%p,K 
% S=zeros(K,N);
% for i=1:N
%     
%     S(:,i)=OMP(x(:,i),A,0.00001);
% end
% 
% %Learn the dictionary using K-SVD algorithm
% for k=1:K
%    
%     %Get non zero sparse codes for each column of dictionary
%     indices=find(S(k,:)~=0);
%     x_iter=x(:,indices);
%     S_iter=S(:,indices);
%     
%     E=x_iter-A*S_iter+A(:,k)*S_iter(k,:);
%     [U_svd,S_svd,V_svd]=svd(E);
%     A(:,k)=U_svd(:,1);
%     S_iter(k,:)=S_svd(1,1)*V_svd(:,1)';
%     S(:,indices)=S_iter;
%     
% end


% D=A;
% m=10; %range of values of m=[10,20,30,40,50,70,90]
% mu=0;
% f=0.001; %range of values of m=[0.001,0.01,0.02,0.05,0.1,0.3]
% sigma=f*sum(sum(abs(phi*x),1))/(m*N);
% 
% 
% phi=(randi([0,1],m,p)-0.5)/(0.5*sqrt(m));
% noise=normrnd(mu,sigma,m,N);
% y=phi*x+noise;
% A=phi*D;
% S=zeros(K,N);
% 
% for i=1:N
%     disp(i);
%     S(:,i)=OMP(y(:,i),A,sigma);
% end
% x_pred=D*S;
% avg_rel_error=0;
% for i=1:N
%    avg_rel_error=avg_rel_error+(norm(x(:,i)-x_pred(:,i))/norm(x(:,i)));
%    
% end
% avg_rel_error=avg_rel_error/N;
% 
% 
% 
% %% PART B
% 
% 
