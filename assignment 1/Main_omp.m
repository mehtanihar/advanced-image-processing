clc;    
clear all;
close all;  
workspace;

nFrames=3;
filename='cars.avi';
noise=1;

video = mmread(filename,1:nFrames);
% height=video.height;
% width=video.width;
height=120;
width=240;

images=zeros(height,width,nFrames);

%% Coded Snapshot

D = dctmtx(64);
coded_pattern=randi([0 1],height,width,nFrames);

for k =1:nFrames
    Image = video.frames(k).cdata;
    I=rgb2gray(Image);
    I=I(end-height+1:end,end-width+1:end);
    images(:,:,k)= double(I);
end

coded_snapshot=zeros(height,width);
for k=1:nFrames
   coded_snapshot=coded_snapshot+coded_pattern(:,:,k).*images(:,:,k); 
end

coded_snapshot=coded_snapshot/255;

if noise==1
    coded_snapshot=imnoise(coded_snapshot,'gaussian',0,4);
end
imshow(coded_snapshot);
    
%% DCT Matrix

dct_matrix=zeros(64*nFrames,64*nFrames);
for k=1:nFrames
   dct_matrix(((k-1)*64+1:(k-1)*64+64),((k-1)*64+1:(k-1)*64+64))=D; 
    
end
disp('DCT');
disp(size(dct_matrix));

%% Non-Overlapping patches------------------------------------
final_images=images;


for i=1:(height/8)
        disp(i)
    for j=1:(width/8)
        C=[];
        
        %Patch codes 
        patch_pattern=coded_pattern((8*(i-1)+1:8*(i-1)+8),(8*(j-1)+1:8*(j-1)+8),:);
        for k=1:nFrames
            a=patch_pattern(:,:,k)';
            C=[C diag(a(:))];
        end
                
        A=C*dct_matrix; 
        
        %Coded snapshot patches
        coded_patch=coded_snapshot((8*(i-1)+1:8*(i-1)+8),(8*(j-1)+1:8*(j-1)+8));
        coded_patch=coded_patch';
        y=coded_patch(:);
        
        %Get sparsest theta using OMP algorithm
        [ theta ] = OMP(y,A);
        
        %Get back the image
        img=dct_matrix*theta;
        
        for k=1:nFrames
            final_images((8*(i-1)+1:8*(i-1)+8),(8*(j-1)+1:8*(j-1)+8),k)=reshape(img((k-1)*64+1:(k-1)*64+64),[8,8])';        
        end
        
    end
end
%% playing the reconstructed movie
cmap = gray(256);
for u=1:nFrames
    % convert the image to a frame
    % mov(u) = im2frame(images{u},cmap);
    mov(u)=im2frame(uint8(final_images(:,:,u)*255), cmap);
    
end
movie(mov);

%% Overlapping patches--------------------------------------------

new_images=zeros(height,width,nFrames);

for i=1:height-7
        disp(i)
    for j=1:width-7
        C=[];
        %Patch codes
        patch_pattern=coded_pattern(i:i+7,j:j+7,:);
        for k=1:nFrames
            a=patch_pattern(:,:,k)';
            C=[C diag(a(:))];
        end
        
        A=C*dct_matrix; 
        
        %Coded snapshot patches
        coded_patch=coded_snapshot(i:i+7,j:j+7);
        coded_patch=coded_patch';
        y=coded_patch(:);
        
        %Get sparsest theta using OMP algorithm
        [ theta ] = OMP(y,A);
        
        %Get back the image
        img=dct_matrix*theta;
        
        for k=1:nFrames
            new_images(i:i+7,j:j+7,k)=new_images(i:i+7,j:j+7,k)+reshape(img((k-1)*64+1:(k-1)*64+64),[8,8])';                
        end
        
    end
end

%Count of overlapping pixels
overlapping_counts=[1 2 3 4 5 6 7 ones(1,height-14)*8 7 6 5 4 3 2 1]'*[1 2 3 4 5 6 7 ones(1,weight-14)*8 7 6 5 4 3 2 1];

%Averaging of pixels
for k=1:nFrames
   new_images(:,:,k)=new_images(:,:,k)./overlapping_counts; 
end

imshow(new_images(:,:,1)/max(max(new_images(:,:,1))));
%% playing the reconstructed movie
cmap = gray(256);
 for u=1:length(images)
     % convert the image to a frame
     % mov(u) = im2frame(images{u},cmap);
       mov(u)=immovie(images{u}, cmap);

 end
%movie(mov,1);




















