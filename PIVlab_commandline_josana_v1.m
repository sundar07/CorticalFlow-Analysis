function PIVlab_commandline_josana_v1(specifyframes)

% Example script how to use PIVlab from the commandline
% You can adjust the settings in "s" and "p", specify a mask and a region of interest

Folder_path = cd;

%Create save directory
if exist(strcat(Folder_path,'/Analysis_v1/Results_PIVlab'),'dir') ~= 7 %check if the directory exists
   mkdir(strcat(Folder_path,'/Analysis_v1/Results_PIVlab'));
end


% threshold the first frame
I = imread(strcat(Folder_path,'/Parameters/A0.tif'));
% This gaussian just blurs the whole image
h=fspecial('gaussian', 10, 5);
% h=fspecial('gaussian', 80, 40);
Ifilt = imfilter(I,h);
It = im2bw(Ifilt,graythresh(Ifilt));

% Subtraction of the thresholded image from the dilated one, gives the outline of the embryo approximately
% Not perfect though
I1 = (imdilate(It,strel('disk',2)) - It);
If = imfill(I1,'holes');
% since the image has been slightly dilated, it is necessary to erode
% once to partially reconstruct the thresholded image
I2 = imerode(If,strel('disk',2));

% Remove anything with area less than 40 pixels. this basically removes
% the background around the embryo
bw = bwareaopen(I2,40);

Ilabel = bwlabel(bw);
Iarea = regionprops(Ilabel, 'area');
D = [Iarea.Area];
% find biggest object and assign as the embryo
[~,ind] = max(D);
%gives the index of where the embryo's located
[i] = find(Ilabel==ind);  
%using this index create a mask that can be applied to the original image
PIVlab_res.mask = nan(size(I));
PIVlab_res.mask(i) = 1;
mask = zeros(size(I));
mask(i) = 1;
mask_peri = bwperim(mask);
PIVlab_res.I_threshold_overlay = imoverlay(I,mask_peri);

close all
figure,imshow(PIVlab_res.I_threshold_overlay)
saveas(gcf,strcat(Folder_path,'/Analysis_v1/Results_PIVlab/Embryo_border'),'tif')

[m n] = find(PIVlab_res.mask == 1);
PIVlab_res.right = max(max(n));
PIVlab_res.left = min(min(n));
PIVlab_res.top = min(min(m));
PIVlab_res.bottom = max(max(m));

PIVlab_res.mfile = 'PIVlab_commandline_paper_v1';

%% Create list of images inside specified directory
dirOutput = dir(strcat(Folder_path,'/*.tif'));
fileNames = {dirOutput.name}';

% Create image sequence array
frame_inc = 1;
for p = specifyframes
    sequence(:,:,frame_inc) = double(imread(fileNames{1},p)) .* PIVlab_res.mask; 
    frame_inc = frame_inc+1;
end

%% Standard PIV Settings
s = cell(10,2); % To make it more readable, let's create a "settings table"
%Parameter                       %Setting           %Options
s{1,1}= 'Int. area 1';           s{1,2}=64;         % window size of first pass
s{2,1}= 'Step size 1';           s{2,2}=32;         % step of first pass
s{3,1}= 'Subpix. finder';        s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
s{4,1}= 'Mask';                  s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
s{5,1}= 'ROI';                   s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
s{6,1}= 'Nr. of passes';         s{6,2}=2;          % 1-4 nr. of passes
s{7,1}= 'Int. area 2';           s{7,2}=32;         % second pass window size
s{8,1}= 'Int. area 3';           s{8,2}=16;         % third pass window size
s{9,1}= 'Int. area 4';           s{9,2}=16;         % fourth pass window size
s{10,1}='Window deformation';    s{10,2}='*linear'; % '*spline' is more accurate, but slower

%% PIV analysis loop

x = cell(length(specifyframes)-1,1);
y = x;
u = x;
v = x;
typevector = x;

% Settings
umin = -12; % minimum allowed u velocity
umax = 12; % maximum allowed u velocity
vmin = -12; % minimum allowed v velocity
vmax = 12; % maximum allowed v velocity
stdthresh = 6; % threshold for standard deviation check
epsilon = 0.15; % epsilon for normalized median test
thresh = 3; % threshold for normalized median test

counter = 0;
PIVlab_res.u_filt = cell(length(specifyframes)-1,1);
PIVlab_res.v_filt = PIVlab_res.u_filt;
PIVlab_res.typevector_filt = PIVlab_res.u_filt;
    
%% PIV analysis loop:
for im=1:length(specifyframes)-1
    counter = counter+1;
    [x{counter} y{counter} u{counter} v{counter} typevector{counter}] = piv_FFTmulti(sequence(:,:,im),sequence(:,:,im+1),s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2});
    clc
    
    u_filtered = u{im,1};
    v_filtered = v{im,1};
    typevector_filtered = typevector{im,1};
    %vellimit check
    u_filtered(u_filtered<umin) = NaN;
    u_filtered(u_filtered>umax) = NaN;
    v_filtered(v_filtered<vmin) = NaN;
    v_filtered(v_filtered>vmax) = NaN;
    % stddev check
    meanu = nanmean(nanmean(u_filtered));
    meanv = nanmean(nanmean(v_filtered));
    std2u = nanstd(reshape(u_filtered,size(u_filtered,1)*size(u_filtered,2),1));
    std2v = nanstd(reshape(v_filtered,size(v_filtered,1)*size(v_filtered,2),1));
    minvalu = meanu-stdthresh*std2u;
    maxvalu = meanu+stdthresh*std2u;
    minvalv = meanv-stdthresh*std2v;
    maxvalv = meanv+stdthresh*std2v;
    u_filtered(u_filtered<minvalu) = NaN;
    u_filtered(u_filtered>maxvalu) = NaN;
    v_filtered(v_filtered<minvalv) = NaN;
    v_filtered(v_filtered>maxvalv) = NaN;
    % normalized median check
    %Westerweel & Scarano (2005): Universal Outlier detection for PIV data
    [J,I] = size(u_filtered);
    medianres = zeros(J,I);
    normfluct = zeros(J,I,2);
    b=1;
    for c=1:2
        if c==1; velcomp=u_filtered;else;velcomp=v_filtered;end %#ok<*NOSEM>
        for i=1+b:I-b
            for j=1+b:J-b
                neigh=velcomp(j-b:j+b,i-b:i+b);
                neighcol=neigh(:);
                neighcol2=[neighcol(1:(2*b+1)*b+b);neighcol((2*b+1)*b+b+2:end)];
                med=median(neighcol2);
                fluct=velcomp(j,i)-med;
                res=neighcol2-med;
                medianres=median(abs(res));
                normfluct(j,i,c)=abs(fluct/(medianres+epsilon));
            end
        end
    end
    info1=(sqrt(normfluct(:,:,1).^2+normfluct(:,:,2).^2)>thresh);
    u_filtered(info1==1)=NaN;
    v_filtered(info1==1)=NaN;

    typevector_filtered(isnan(u_filtered))=2;
    typevector_filtered(isnan(v_filtered))=2;
    typevector_filtered(typevector{im,1}==0)=0; %restores typevector for mask
    
    %Interpolate missing data
    u_filtered=inpaint_nans(u_filtered,4);
    v_filtered=inpaint_nans(v_filtered,4);
    
    PIVlab_res.u_filt{im,1}=u_filtered;
    PIVlab_res.v_filt{im,1}=v_filtered;
    PIVlab_res.typevector_filt{im,1}=typevector_filtered;

    % Graphical output (disable to improve speed)
    clf;imshow(sequence(:,:,im+1),[])
    hold on
    quiver(x{counter},y{counter},PIVlab_res.u_filt{counter},PIVlab_res.v_filt{counter},0,'g','LineWidth',0.3);
%     title(['Frame ' num2str(im)]);
    export_fig(gcf,strcat(Folder_path,'/Analysis_v1/Results_PIVlab/PIV_frame',num2str(im+1),'.pdf'),'-q101','-m2','-r200')
end

PIVlab_res.x = x; PIVlab_res.y = y;
PIVlab_res.u = u; PIVlab_res.v = v;
save(strcat(Folder_path,'/Analysis_v1/Results_PIVlab/PIVlab_res.mat'),'PIVlab_res')
