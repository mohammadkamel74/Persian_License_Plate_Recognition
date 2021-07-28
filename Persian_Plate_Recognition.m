tic
%load imgfildata;
% The steps are as follows:
% 1. Vertical edge detection
% 2. Histogram analysis
% 3. Vertical and Horizontal Dilation
% 4. Finding Regions in Common
% 5. Horizontal Dilation
% 6. Erosion
% 7. Post processing

%% Read Image
In    = 1234;    %Name of the image file
Im    = imread(fullfile( strcat(num2str(In), '.jpg'))); % read image
tform = affine2d([1 .4 0; 0 1 0; 0 0 1])
Im = imwarp(Im,tform);
I     = im2double(rgb2gray(Im));        % rgb to gray
figure();imshow(I);title('RGB Image')

%% Sobel Masking 
SM    = [-1 0 1;-2 0 2;-1 0 1];         % Sobel Vertical Mask
IS    = imfilter(I,SM,'replicate');     % Filter Image Using Sobel Mask
IS    = IS.^2;                          % Consider Just Value of Edges & Fray Weak Edges
%% Normalization
IS    = (IS-min(IS(:)))/(max(IS(:))-min(IS(:))); % Normalization
%% Threshold (Otsu)
level = graythresh(IS);                 % Threshold Based on Otsu Method
IS    = im2bw(IS,level);
%% Histogram
S     = sum(IS,2);                      % Edge Horizontal Histogram
%% Plate Location
T1    = 0.35;                           % Threshold On Edge Histogram
PR    = find(S > (T1*max(S)));          % Candidate Plate Rows
%% Masked Plate
Msk   = zeros(size(I));
Msk(PR,:) = 1;                          % Mask
MB    = Msk.*IS;                        % Candidate Plate (Edge Image)
figure();imshow(MB);title('4')
%% Morphology (Dilation - Vertical)
Dy    = strel('rectangle',[80,4]);      % Vertical Extension
MBy   = imdilate(MB,Dy);                % By Dilation
MBy   = imfill(MBy,'holes');            % Fill Holes
%% Morphology (Dilation - Horizontal)
Dx    = strel('rectangle',[4,80]);      % Horizontal Extension
MBx   = imdilate(MB,Dx);                % By Dilation
MBx   = imfill(MBx,'holes');            % Fill Holes
%% Joint Places
BIM   = MBx.*MBy;                       % Joint Places
%% Morphology (Dilation - Horizontal)
Dy    = strel('rectangle',[4,30]);      % Horizontal Extension
MM    = imdilate(BIM,Dy);               % By Dilation
MM    = imfill(MM,'holes');             % Fill Holes
%% Erosion
Dr    = strel('line',50,0);             % Erosion
BL    = imerode(MM,Dr);
%% Find Biggest Binary Region (As a Plate Place)
[L,num] = bwlabel(BL);                  % Label (Binary Regions)               
Areas   = zeros(num,1);
for i = 1:num                           % Compute Area Of Every Region
[r,c,v]  = find(L == i);                % Find Indexes
Areas(i) = sum(v);                      % Compute Area    
end 
try
    [La,Lb] = find(Areas==max(Areas));      % Biggest Binary Region Index
    %% Post Processing
    [a,b]   = find(L==La);                  % Find Biggest Binary Region (Plate)
    [nRow,nCol] = size(I);
    FM      = zeros(nRow,nCol);             % Smooth and Enlarge Plate Place
    T       = 10;                           % Extend Plate Region By T Pixel
    jr      = (min(a)-T :max(a)+T);
    jc      = (min(b)-T :max(b)+T);
    jr      = jr(jr >= 1 & jr <= nRow);
    jc      = jc(jc >= 1 & jc <= nCol);
    FM(jr,jc) = 1; 
    PL      = FM.*I;                        % Detected Plate
    crop = imcrop(Im,[min(jc)-20 min(jr) max(jc)-min(jc)+40 max(jr)-min(jr)]);
    Im=imresize(crop,2);
    I     = im2double(rgb2gray(Im));        % rgb to gray
     figure();imshow(I);title('5')
    %% Sobel Masking 
    SM    = [-1 0 1;-2 0 2;-1 0 1];         % Sobel Vertical Mask
    IS    = imfilter(I,SM,'replicate');     % Filter Image Using Sobel Mask
    IS    = IS.^2;                          % Consider Just Value of Edges & Fray Weak Edges
    %% Normalization
    IS    = (IS-min(IS(:)))/(max(IS(:))-min(IS(:))); % Normalization
    %% Threshold (Otsu)
    level = graythresh(IS);                 % Threshold Based on Otsu Method
    IS    = im2bw(IS,level);
    %% Histogram
    S     = sum(IS,2);                      % Edge Horizontal Histogram
    %% Plate Location
    T1    = 0.35;                           % Threshold On Edge Histogram
    PR    = find(S > (T1*max(S)));          % Candidate Plate Rows
    %% Masked Plate
    Msk   = zeros(size(I));
    Msk(PR,:) = 1;                          % Mask
    MB    = Msk.*IS;                        % Candidate Plate (Edge Image)
    figure();imshow(MB);title('6')
    %% Morphology (Dilation - Vertical)
    Dy    = strel('rectangle',[80,4]);      % Vertical Extension
    MBy   = imdilate(MB,Dy);                % By Dilation
    MBy   = imfill(MBy,'holes');            % Fill Holes
    %% Morphology (Dilation - Horizontal)
    Dx    = strel('rectangle',[4,80]);      % Horizontal Extension
    MBx   = imdilate(MB,Dx);                % By Dilation
    MBx   = imfill(MBx,'holes');            % Fill Holes
    %% Joint Places
    BIM   = MBx.*MBy;                       % Joint Places
    %% Morphology (Dilation - Horizontal)
    Dy    = strel('rectangle',[4,30]);      % Horizontal Extension
    MM    = imdilate(BIM,Dy);               % By Dilation
    MM    = imfill(MM,'holes');             % Fill Holes
    %% Erosion
    Dr    = strel('line',50,0);             % Erosion
    BL    = imerode(MM,Dr);
    %% Find Biggest Binary Region (As a Plate Place)
    [L,num] = bwlabel(BL);                  % Label (Binary Regions)               
    Areas   = zeros(num,1);
    for i = 1:num                           % Compute Area Of Every Region
    [r,c,v]  = find(L == i);                % Find Indexes
    Areas(i) = sum(v);                      % Compute Area    
    end 
    [La,Lb] = find(Areas==max(Areas));      % Biggest Binary Region Index
    %% Post Processing
    [a,b]   = find(L==La);                  % Find Biggest Binary Region (Plate)
    [nRow,nCol] = size(I);
    FM      = zeros(nRow,nCol);             % Smooth and Enlarge Plate Place
    T       = 10;                           % Extend Plate Region By T Pixel
    jr      = (min(a)-T :max(a)+T);
    jc      = (min(b)-T :max(b)+T);
    jr      = jr(jr >= 1 & jr <= nRow);
    jc      = jc(jc >= 1 & jc <= nCol);
    FM(jr,jc) = 1; 
    PL      = FM.*I;                        % Detected Plate
    %% Plot
    crop = imcrop(Im,[min(jc) min(jr) max(jc)-min(jc) max(jr)-min(jr)]);
    if max(jc)-min(jc)>200
        plate_status=0
    elseif max(jc)-min(jc)<100
        plate_status=0
    elseif max(jr)-min(jr)>70
        plate_status=0
    elseif max(jr)-min(jr)<30
        plate_status=0
    else
        plate_status=1
    end
catch exception
   plate_status=0
end

if plate_status==1
    %% Image segmentation and extraction
    A=rgb2gray(crop);
    B=imresize(A,[100,300]);
    pic=B;
    C = imadjust(B);
    threshold = graythresh(C);
    imagen =~imbinarize(C,(threshold));
    %% Remove all object containing fewer than 200 pixels
    imagen = bwareaopen(imagen,200);
    %% Remove all object containing more than 1200 pixels
    cc = bwconncomp(imagen);
    stats = regionprops(cc);
    limit = 1200;
    removeMask = [stats.Area]>limit;
    imagen(cat(1,cc.PixelIdxList{removeMask})) = false;
    %% Remove extra line
    var=0;                                                   %#7
    for i=1:size(imagen,1)
        for j=1:size(imagen,2)
            if (imagen(i,j)==1 && var==0)
               index1i=i;
               var=1;
            end
            if (imagen(i,j)==1 )
               index2i=i;     
            end
        end
    end
    figure;imshow(imagen);title('GrayScalePlate');
    figure;imshow(pic);title('ROIPlate');
        J = imcrop(B,[1 index1i  size(imagen,2) index2i-index1i ]);
        %#8
        C = imadjust(J);
        threshold = graythresh(C);
        imagen =~imbinarize(C,(threshold));
        %% Remove all object containing fewer than 200 pixels
        imagen = bwareaopen(imagen,200);
        %% Remove all object containing more than 1200 pixels
        cc = bwconncomp(imagen);
        stats = regionprops(cc);
        limit2 = 1200;
        removeMask = [stats.Area]>limit2;
        imagen(cat(1,cc.PixelIdxList{removeMask})) = false;
        
        relese2=0;                                        %#3
        r=100;
        c=100;
        KK= bwlabel(imagen);
        for i=1:10
             for j=1:10
                 if KK(i,j)==1
                     relese2=1;
                 end
             end
        end
        if relese2==1
               [r, c] = find(KK==1);
               rc = [r c];
        end
        for i=1:size(r,1)
            imagen(r(i,1),c(i,1))=0;
        end
        %#4
        relese3=0;                                        %#5
        r=1;
        c=1;
        KK= bwlabel(imagen);
        for i=1:20
             for j=size(imagen,2)-5:size(imagen,2)
                 if KK(i,j)==max(KK(:))
                     relese3=1;
                 end
             end
        end
        if relese3==1
               [r, c] = find(KK==max(KK(:)));
               rc = [r c];
        end
        for i=1:size(r,1)
            imagen(r(i,1),c(i,1))=0;
        end
        %#6
        %% Label connected components
        [L Ne]=bwlabel(imagen);
        %% Measure properties of image regions
        propied=regionprops(L,'BoundingBox');
        %hold on
        %% Objects extraction and remove 2 unconnected part
        %figure
        for n=1:Ne                                          %#1
          relese=0;
          [r,c] = find(L==n);
          n1=imagen(min(r):max(r),min(c):max(c));
          LL= bwlabel(n1);
          for i=1:size(n1,1)
             for j=1:size(n1,2) 
                 if LL(i,j)==2
                     relese=1;
                 end
             end
          end
          if relese==0
            %imshow(~n1);
            imwrite(n1, sprintf('file%.02d.jpg', n));
            %pause(0.5)                                      %#2
          end
        end
end
final_output=[];
t=[];
srcFiles = dir('E:\UNIVERSITY\PROJECTS\PORREZA_LAB\MACHINE\plate-whole-project\complete\*.jpg');
for n=1:length(srcFiles)
  filename = strcat('E:\UNIVERSITY\PROJECTS\PORREZA_LAB\MACHINE\plate-whole-project\complete\',srcFiles(n).name);  
  n1=imread(filename);
  n1=imresize(n1,[42,24]);
  x=[ ];
totalLetters=size(imgfile,2);
 for k=1:totalLetters
    y=corr2(imgfile{1,k},n1);
    x=[x y];
 end
 t=[t max(x)];
 if max(x)>.45
 z=find(x==max(x));
 out=cell2mat(imgfile(2,z));
final_output=[final_output out];
end
end
file = fopen('number_Plate.txt', 'wt');
    fprintf(file,'%s\n',final_output);
    fclose(file);                     
toc
clc