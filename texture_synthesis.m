for or=4:4
%read in image
order=mat2str(or);
I=im2double(imread(['./texture/texture',order,'.jpg']));
sizeI_1=size(I,1);
sizeI_2=size(I,2);

%threshold
window_half=ceil(size(I,1)/3);
L_size=ceil(size(I,1)/2)*3;%大图size
sigma=window_half/3;



%devide RGB and grey
if size(size(I),2)==3
    RGB=3;
else
    RGB=1;
    I(:,:,2)=zeros(size(I,1),size(I,2));
    I(:,:,3)=zeros(size(I,1),size(I,2));
end


%sample convolution by window_half*2+1
I_convol=[];
for i=1:RGB
I_convol(:,:,i) = im2col(I(:,:,i), [window_half*2+1 window_half*2+1],'sliding');
end


%Generate Large, Left up coner:sample
Large=[];
for i=1:3
Large(:,:,i)=-ones(L_size,L_size);
end
Large(1:sizeI_1, 1:sizeI_2,:) = I;


%add window_half+1 size -2 around Large and filled for easy culculation
Large=[-2*ones(window_half+1,L_size,3); Large];
Large=[Large; -2*ones(window_half+1,L_size,3)];
Large=[Large -2*ones(size(Large,1),window_half+1,3)];
Large=[-2*ones(size(Large,1),window_half+1,3) Large];


%Generate Gaussian parameter(initial)
Gaus_ini=fspecial('gaussian',window_half*2+1, sigma);
Gaus_ini=Gaus_ini(:);


%Fill each layer around the sample in turn
if sizeI_1>=sizeI_2
    size_min=sizeI_2;
else
    size_min=sizeI_1;
end

for ftimes=window_half+2+size_min:size(Large,1)-1-window_half
        %get pixels need to be filles
        fill_list=[];
        for i=window_half+2:ftimes-1
            fill_list=[fill_list;i ftimes];
        end
        for j=window_half+2:ftimes-1
            fill_list=[fill_list;ftimes j];
        end
        fill_list=[fill_list;ftimes ftimes];
        for n=1:size(fill_list,1)
            x=fill_list(n,1);
            y=fill_list(n,2);
            [window_xy,Gaus]=window_gene(x,y,Large,Gaus_ini,window_half,RGB);
            matchpoint=match(window_xy,Gaus,I_convol,RGB);
            for d=1:RGB
            Large(x,y,d)=matchpoint(d,1);
            end
        end
       
        
    
end


%write out image
L_s=size(Large,1)+1;
    for i=1:window_half+1
         Large(1,:,:)=[];
         Large(:,1,:)=[];
         Large(L_s-i*2,:,:)=[];
         Large(:,L_s-i*2,:)=[];  
    end
    if RGB==1
        Largenew=[];
        Largenew=Large(:,:,1);
        Large=[];
        Large=Largenew;
    end
Large=im2uint8(Large);
imwrite(Large,['new',order,'.png']);
end


%generate window vector aroung x,y, and modify Gaus_ini by -2 and -1
function [window_xy,Gaus]=window_gene(x,y,Large,Gaus_ini,window_half,RGB)
window_xy=zeros(window_half*2+1,window_half*2+1,3);
Gaus=reshape(Gaus_ini,window_half*2+1,window_half*2+1);
for i=-window_half:window_half
    for j=-window_half:window_half
        for d=1:RGB
        window_xy(i+window_half+1,j+window_half+1,d)=Large(i+x,j+y,d);
        if Large(i+x,j+y,d)==-1 || Large(i+x,j+y,d)==-2
            Gaus(i+window_half+1,j+window_half+1)=0;
        end
        end
    end
end
window_1=window_xy(:,:,1);
window_2=window_xy(:,:,2);
window_3=window_xy(:,:,3);
window_xy=[];
window_xy(:,:,1)=window_1(:);
window_xy(:,:,2)=window_2(:);
window_xy(:,:,3)=window_3(:);

Gaus=Gaus(:);
end


%find match for x,y accoring to window_xy
function matchpoint=match(window_xy,Gaus,I_convol,RGB)
min=inf;
matchpoint_ca=[];
matchpoint=zeros(3,1);
for i=1:size(I_convol,2)
    sum=0;
    for d=1:RGB
    sum= sum+((window_xy(:,:,d)-I_convol(:,i,d)).*(window_xy(:,:,d)-I_convol(:,i,d)))'*Gaus;
    end
    if sum<min
        matchpoint=[];
        min=sum;
        matchpoint_ca=[ceil(size(I_convol,1)/2) i];
    elseif sum==min
        matchpoint_ca=[matchpoint_ca;ceil(size(I_convol,1)/2) i];
    end
    A=randsample(size(matchpoint_ca,1),1,false);
    s=matchpoint_ca(A,:);
    si=s(1,1);
    sj=s(1,2);
    for d=1:RGB
        matchpoint(d,1)=I_convol(si,sj,d);
    end
end
end

