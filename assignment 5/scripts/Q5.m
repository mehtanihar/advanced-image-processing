c=1;
alpha=3;
eigen=[];
i=1:128;
val=c*i.^(-alpha);
D=diag(val);
N=128;
U=rand(128,128);
tolerance=1e-6;
ui = randn(N,1);
U(:,1) = ui ./ norm(ui);  
for i=2:N
  norm_factor = 0;
  while norm_factor<tolerance
    ui = randn(N,1);
    ui = ui-U(:,1:i-1)*(U(:,1:i-1).'*ui);
    norm_factor = norm(ui);
  end
  U(:,i) = ui ./ norm_factor;
end 

sigma=U*D*U';
x=mvnrnd(zeros(N,1),sigma,10)';
m_vec=[40,50,64,80,100,120];
avg_rmse=zeros(6,1);
for i=1:length(m_vec)
    m=m_vec(i);

    phi=normrnd(0,sqrt(1/m),[m,N]);
    sigma_small=0.01*mean(mean(abs(phi*x)));
    noise=mvnrnd(zeros(m,1),sigma_small*eye(m),10)';
    y=phi*x+noise;

    x_recons=(phi'*phi+sigma_small^2*inv(sigma))\phi'*y;
    avg_rmse(i,1)=mean(sqrt(mean((x-x_recons).^2,1)));

end
plot(m_vec,avg_rmse);




