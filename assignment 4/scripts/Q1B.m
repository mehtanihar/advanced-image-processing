images=loadMNISTImages('../data/train-images.idx3-ubyte');
indices=randi([1 60000],[600 1]);
images=images(:,indices);
N=600;
p=16*16;
m=16; %range of values of m=[10,20,30,40,50,70,90]
mu=0;
f=0.01; %range of values of m=[0.001,0.01,0.02,0.05,0.1,0.3]
sigma=0;
K=20;
D=randn(p,K);
%Unit normalize columns of D
for i=1:K
   D(:,i)=D(:,i)/norm(D(:,i));
end

x=zeros(p,N);
for i=1:N
    image=reshape(images(:,i),[28,28]);
    image_new=bilinearInterpolation(image,[16,16]);
    x(:,i)=reshape(image_new,[p,1]);
end

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
        S(:,i)=OMP(y(:,i),reshape(A(i,:,:),[m,K]),sigma);
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

        dk=G\h;
        
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
