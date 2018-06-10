%% PART A: FILTERED BACK PROJECTION

image1=double(imread('../input/slice1.png'));
[row,col,dim]=size(image1);
angles=round(rand(18,1)*180);
radon_transform=radon(image1,angles);
FBP_reconstructed_image1=iradon(radon_transform,angles,'linear','Ram-Lak',1,row);

image2=double(imread('../input/slice2.png'));
[row,col,dim]=size(image2);
angles=round(rand(18,1)*180);
radon_transform=radon(image2,angles);
FBP_reconstructed_image2=iradon(radon_transform,angles,'linear','Ram-Lak',1,row);

figure; subplot(2,2,1); imshow(image1/255); title('Original image 1');
subplot(2,2,2); imshow(FBP_reconstructed_image1/255); title('FBP Reconstructed image 1');
subplot(2,2,3); imshow(image2/255); title('Original image 2');
subplot(2,2,4); imshow(FBP_reconstructed_image2/255); title('FBP Reconstructed image 2');

disp('Error:');
disp(norm(FBP_reconstructed_image1/255-image1/255));
disp(norm(FBP_reconstructed_image2/255-image2/255));

%% PART B: COMPRESSED SENSING
lambda=10;
image1=double(imread('../input/slice1.png'));
[row,col,dim]=size(image1);
radon_transform=radon(image1,angles);
y=radon_transform(:);
A=radon_operator(row,col,angles);
A_T=A';
val=l1_ls(A,A_T,size(y,1),row*col,y,lambda,0.1,'false');
recons_image1b=dct2(reshape(val,[row,col]))/255;

figure; subplot(2,2,1); imshow(image1/255); title('Original image 1');
subplot(2,2,2); imshow(recons_image1b); title('Reconstructed image 1');
disp('Error:');
disp(norm(recons_image1b-image1/255));

lambda=10;
image2=double(imread('../input/slice2.png'));
[row,col,dim]=size(image2);
radon_transform=radon(image2,angles);
y=radon_transform(:);
A=radon_operator(row,col,angles);
A_T=A';
val=l1_ls(A,A_T,size(y,1),row*col,y,lambda,1,'false');
recons_image2b=dct2(reshape(val,[row,col]))/255; %5.8

subplot(2,2,3); imshow(image2/255); title('Original image 2');
subplot(2,2,4); imshow(recons_image2b); title('Reconstructed image 2');
disp('Error:');
disp(norm(recons_image2b-image2/255));

%% PART C: COUPLED COMPRESSED SENSING WITH 2 SLICES

lambda=300;
image1=double(imread('../input/slice1.png'));
arr=1:180;
angles1=round(rand(18,1)*180);
arr(angles1)=[];
[row,col,dim]=size(image1);
radon_transform=radon(image1,angles1);
y1=radon_transform(:);

image2=double(imread('../input/slice2.png'));
angles2=round(rand(18,1)*162);
angles2=arr(angles2)';
radon_transform=radon(image2,angles2);
y2=radon_transform(:);

A=radon_operator2(row,col,angles1,angles2);
A_T=A';
y=[y1;y2];
val=l1_ls(A,A_T,size(y,1),2*row*col,y,lambda,1,'false');
recons_image1c=dct2(reshape(val(1:row*col),[row,col]))/255;
recons_image2c=recons_image1c+dct2(reshape(val(row*col+1:end),[row,col]))/255; %6.3

figure; subplot(2,2,1); imshow(image1/255); title('Original image 1');
subplot(2,2,2); imshow(recons_image1c); title('Reconstructed image 1');
subplot(2,2,3); imshow(image2/255); title('Original image 2');
subplot(2,2,4); imshow(recons_image2c); title('Reconstructed image 2');

disp('Error in image 1:');
disp(norm(recons_image1c-image1/255));
disp('Error in image 2:');
disp(norm(recons_image2c-image2/255));

%% PART D: COUPLED COMPRESSED SENSING WITH 3 SLICES

lambda=500;
pad=0;
image1_old=double(imread('../input/slice_50.png'));
image2_old=double(imread('../input/slice_51.png'));
image3_old=double(imread('../input/slice_52.png'));

[row,col,dim]=size(image1_old);
if(row<col)
    pad=round((col-row)/2);
    image1=zeros(2*pad+row,col);
    image2=zeros(2*pad+row,col);
    image3=zeros(2*pad+row,col);
    image1(pad+1:pad+row,1:col)=image1_old;
    image2(pad+1:pad+row,1:col)=image2_old;
    image3(pad+1:pad+row,1:col)=image3_old;
    
elseif(row>col)
    pad=round((row-col)/2);
    image1=zeros(row,2*pad+col);
    image2=zeros(2*pad+row,2*pad+col);
    image3=zeros(row,2*pad+col);
    image1(1:row,pad+1:pad+col)=image1_old;
    image2(1:row,pad+1:pad+col)=image2_old;
    image3(1:row,pad+1:pad+col)=image3_old;
else
    image1=image1_old;
    image2=image2_old;
    image3=image3_old;
end

[row,col,dim]=size(image1);

arr=1:180;
angles1=round(rand(18,1)*180);
arr(angles1)=[];

radon_transform=radon(image1,angles1);
y1=radon_transform(:);


indices=round(rand(18,1)*162);
angles2=arr(indices)';
arr(indices)=[];
radon_transform=radon(image2,angles2);
y2=radon_transform(:);


indices=round(rand(18,1)*144);
angles3=arr(indices)';
arr(indices)=[];
radon_transform=radon(image3,angles3);
y3=radon_transform(:);

A=radon_operator3(row,col,angles1,angles2,angles3);
A_T=A';
y=[y1;y2;y3];

val=l1_ls(A,A_T,size(y,1),3*row*col,y,lambda,4,'false');

recons_image1d=dct2(reshape(val(1:row*col),[row,col]))/255;
recons_image2d=recons_image1d+dct2(reshape(val(row*col+1:2*row*col),[row,col]))/255;
recons_image3d=recons_image2d+dct2(reshape(val(2*row*col+1:3*row*col),[row,col]))/255;

figure; subplot(3,2,1); imshow(image1/255); title('Original slice 1');
subplot(3,2,2); imshow(recons_image1d); title('Reconstructed slice 1');
subplot(3,2,3); imshow(image2/255); title('Original slice 2');
subplot(3,2,4); imshow(recons_image2d); title('Reconstructed slice 2');
subplot(3,2,5); imshow(image3/255); title('Original slice 3');
subplot(3,2,6); imshow(recons_image3d); title('Reconstructed slice 3');

disp('Error in image 1:');
disp(norm(recons_image1d-image1/255));
disp('Error in image 2:');
disp(norm(recons_image2d-image2/255));
disp('Error in image 3:');
disp(norm(recons_image3d-image3/255));

