rng(0,'twister');
orig_image=imread('../input/barbara256.png');
[H,W]=size(orig_image);
orig_image=double(orig_image);
sigma = 10;
mean = 0;
y = sigma*randn(H,W) + mean;
image=orig_image+y;
imshow(image/255)
image=image;
%% PART A: Denoising


U=kron(dctmtx(8)',dctmtx(8)')';
stride=1;
A=U;
alpha=1.5*eigs(A'*A,1);
lambda=1;
final_image=zeros(H,W);

for i=1:stride:H-7
    disp(i)
    for j=1:stride:W-7
        
        y=image(i:i+7,j:j+7);
        
        y=y';
        y=y(:);
        k=0;
        theta=rand(64,1);
        prev_theta=theta*0;
        count=0;
        
        while(norm(theta-prev_theta)>0.00001)
           prev_theta=theta;
           count=count+1;
           theta = soft(theta + (A'*(y - A*theta))/alpha,lambda/(2*alpha));
           
           
        end
                
        final_image(i:i+7,j:j+7)=final_image(i:i+7,j:j+7)+reshape(A*theta,[8,8])';
       
    end
end

%Average overlapping pixels
avg_row_vector=zeros(1,W);
avg_col_vector=zeros(H,1);
for i=1:stride:H-7
    avg_col_vector(i:i+7,1)=avg_col_vector(i:i+7,1)+ones(8,1);
end

avg_col_vector(i+8:end,1)=1;

for i=1:stride:W-7
    avg_row_vector(1,i:i+7)=avg_row_vector(1,i:i+7)+ones(1,8);
end

avg_row_vector(1,i+8:end)=1;

avg_matrix=avg_col_vector*avg_row_vector;
final_image=final_image./avg_matrix;

imshow(final_image/255);

%% PART B: Image reconstruction


image=orig_image/255;
U=kron(dctmtx(8)',dctmtx(8)')';
phi=rand(32,64);
stride=4;
A=phi*U;
alpha=1.5*eigs(A'*A,1);
lambda=1;
final_image=zeros(H,W);
for i=1:stride:H-7
    disp(i)
    for j=1:stride:W-7
        
        x_true=image(i:i+7,j:j+7);        
        x_true=x_true';
        x_true=x_true(:);
        y=phi*x_true;
       
        theta=rand(64,1);
        prev_theta=theta*0;
        count=0;
        
        while(norm(theta-prev_theta)>0.00001)
           prev_theta=theta;
           count=count+1;
           theta = soft(theta + (A'*(y - A*theta))/alpha,lambda/(2*alpha));
           %disp(norm(theta-prev_theta))
        end
        
        
        
        final_image(i:i+7,j:j+7)=final_image(i:i+7,j:j+7)+reshape(U*theta,[8,8])';
       
        %disp(norm(theta-prev_theta));
        %disp(norm(final_image(i:i+7,j:j+7)-orig_image(i:i+7,j:j+7)/255.0))
    end
end

%Average overlapping pixels
avg_row_vector=zeros(1,W);
avg_col_vector=zeros(H,1);
for i=1:stride:H-7
    avg_col_vector(i:i+7,1)=avg_col_vector(i:i+7,1)+ones(8,1);
end

avg_col_vector(i+8:end,1)=1;

for i=1:stride:W-7
    avg_row_vector(1,i:i+7)=avg_row_vector(1,i:i+7)+ones(1,8);
end

avg_row_vector(1,i+8:end)=1;

avg_matrix=avg_col_vector*avg_row_vector;
final_image=final_image./avg_matrix;

imshow(final_image/255);


%% PART C: Convolved signal

x=zeros(1,100);
x(randi(100,[1,10]))=rand([1,10]);
h=[1 2 3 4 3 2 1]/16;
y=x;
for i=1:100
    disp(i)
    if(i>=4 && i<=97)
        y(i)=h*x(1,i-3:i+3)';
    elseif(i<4)
        y(i)=h(5-i:7)*x(1,1:i+3)';
    elseif(i>97)
        y(i)=h(1:104-i)*x(1,i-3:100)';
    end
    
end

sigma = 0.05*norm(x);
mean = 0;
noise = sigma*randn(1,100) + mean;
y=y+noise;

A=[];
for i=1:100
   if(i<4)
      A=[A; [h(5-i:7) zeros(1,97-i)]]; 
   elseif(i>=4 && i<=97)
       A=[A; [zeros(1,i-4) h zeros(1,97-i)]];
   elseif(i>97)
       A=[A; [zeros(1,i-4) h(1:104-i)]];
   end
    
end


alpha=6*eigs(A'*A,1);
lambda=1;
theta=rand(100,1);

prev_theta=100*theta;
count=0;
while(norm(theta-prev_theta)>0.1)
   
   prev_theta=theta;
   count=count+1;
   theta = soft(theta + (A'*(y' - A*theta))/alpha,lambda/(2*alpha));
   disp(norm(theta));
end

disp(norm(theta'-x));