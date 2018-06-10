b=1;
dx=0.1;
x=dx:dx:b;
image=im2double(imread('../input/d6.png'));
D = dctmtx(8);
[row,col,dim]=size(image);
f=[];
count=0;
for i=1:8:row-7
    for j=1:8:col-7
        disp(i);
        A=image(i:i+7,j:j+7);
        dct_coeff=dct2(A);
        count=count+1;
        f=[f dct_coeff(:)];
    end
end

%% 

len=size(f,2);
p=zeros(64,len);
for index=1:64
    disp(index);
    for i=1:len
        
        p(index,i)=sum(sqrt(2)/(sqrt(pi))*exp(-f(index,i)^2./(2*x.^2))*dx);
        %p(index,i)=1/sqrt(2)*exp(-sqrt(2)*abs(f(index,i)));%Original DCT
    end
end

%% 

hist(p(32,:));


