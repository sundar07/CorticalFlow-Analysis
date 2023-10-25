function Bin_vel_josana_v1(x_rot,y_rot,vx_rot,vy_rot,No_bins,specifyframes,Ppole,Apole,top,bottom,pixelsize,frameinterval,half_stripewidth)

% function to bin the embryo into columns, find the mean velocity of
% these No_bins and plot them

% Get the posterior and anterior pole pixels from Rot_vec_final3 and find
% the embryo length

% Version 4 - abs calculations and bin calculations have been removed. We decided to do analysis only in stripes

Folder_path = cd;

% convert pixels/frame to um/min
vx_rot = vx_rot.*60.*pixelsize./frameinterval;
vy_rot = vy_rot.*60.*pixelsize./frameinterval;

elength = Ppole - Apole;
Bin.ecenter_y = (bottom - top)/2;
binwidth = elength/No_bins; % This gives binwidth of each bin
Bin.No_binseq = round(Apole:binwidth:Ppole);

% No_frames = length(specifyframes)-1;
No_frames = length(specifyframes);

% Save all the parameters
Bin.Parameters.No_bins = No_bins;
Bin.Parameters.specifyframes = specifyframes;
Bin.mfile = 'Bin_vel_josana_v1';

%%
% Pre-allocate the variables
Bin.vx_prof = NaN(No_frames,No_bins);
Bin.vy_prof = NaN(No_frames,No_bins);
Bin.mag_prof = NaN(No_frames,No_bins);
Bin.Remove.vx_stripe = NaN(No_frames,No_bins);
Bin.Remove.vy_stripe = NaN(No_frames,No_bins);
Bin.Remove.mag_prof = NaN(No_frames,No_bins);

[post_mag_chc] = find(x_rot >= Bin.No_binseq(14) & x_rot <= Bin.No_binseq(16) & y_rot >= (top+Bin.ecenter_y-half_stripewidth/2) & y_rot <= (top+Bin.ecenter_y+half_stripewidth/2));
[post] = find(x_rot >= Bin.No_binseq(13) & x_rot <= Bin.No_binseq(17) & y_rot >= (top+Bin.ecenter_y-half_stripewidth) & y_rot <= (top+Bin.ecenter_y+half_stripewidth));
[ant] = find(x_rot >= Bin.No_binseq(3) & x_rot <= Bin.No_binseq(7) & y_rot >= (top+Bin.ecenter_y-half_stripewidth) & y_rot <= (top+Bin.ecenter_y+half_stripewidth));
[mid] = find(x_rot >= Bin.No_binseq(8) & x_rot <= Bin.No_binseq(12) & y_rot >= (top+Bin.ecenter_y-half_stripewidth) & y_rot <= (top+Bin.ecenter_y+half_stripewidth));

% af stands for all frames
Bin.vxp_af = nan(1,No_frames);
Bin.vyp_af = nan(1,No_frames);
Bin.magp_af = nan(1,No_frames);
Bin.magp_chc_af = nan(1,No_frames);

Bin.vxa_af = nan(1,No_frames);
Bin.vya_af = nan(1,No_frames);
Bin.maga_af = nan(1,No_frames);

Bin.vxm_af = nan(1,No_frames);
Bin.vym_af = nan(1,No_frames);
Bin.magm_af = nan(1,No_frames);

%%

for frame = 1:No_frames
    
    vxf = vx_rot(:,:,frame); % choose the first frame
    vyf = vy_rot(:,:,frame);
    [~, magf] = cart2pol(vxf,vyf);
    
    Bin.vxp_af(1,frame) = nanmean(vxf(post));
    Bin.vyp_af(1,frame) = nanmean(vyf(post));
    Bin.magp_af(1,frame) = nanmean(magf(post));
    Bin.magp_chc_af(1,frame) = nanmean(magf(post_mag_chc));
    
    Bin.vxa_af(1,frame) = nanmean(vxf(ant));
    Bin.vya_af(1,frame) = nanmean(vyf(ant));
    Bin.maga_af(1,frame) = nanmean(magf(ant));
    
    Bin.vxm_af(1,frame) = nanmean(vxf(mid));
    Bin.vym_af(1,frame) = nanmean(vyf(mid));
    Bin.magm_af(1,frame) = nanmean(magf(mid));
    
    for bin = 1:No_bins % Run the loop for each bin
        
        [stripe] = find(x_rot > Bin.No_binseq(bin) & x_rot <= Bin.No_binseq(bin+1) & y_rot >= (top+Bin.ecenter_y-half_stripewidth) & y_rot <= (top+Bin.ecenter_y+half_stripewidth));
        
        %%
        % get all the velocities for the positions that you have chosen in the previous step
        vx_stripe = vxf(stripe);
        vy_stripe = vyf(stripe);
        mag_stripe = magf(stripe);
        
        % Find mean and standard deviation and use this to filter noisy velocities
        vx_mean_stripe = nanmean(vx_stripe);
        vy_mean_stripe = nanmean(vy_stripe);
        mag_mean_stripe = nanmean(mag_stripe);
        vx_std_stripe = nanstd(vx_stripe);
        vy_std_stripe = nanstd(vy_stripe);
        mag_std_stripe = nanstd(mag_stripe);
        
        % Remove velocities which are greater than 2 SD
        [vx_inc] = find(vx_stripe >= (vx_mean_stripe-2*vx_std_stripe) & vx_stripe <= (vx_mean_stripe+2*vx_std_stripe));
        Bin.vx_prof_af(frame,bin) = nanmean(vx_stripe(vx_inc));
        [vy_inc] = find(vy_stripe >= (vy_mean_stripe-2*vy_std_stripe) & vy_stripe <= (vy_mean_stripe+2*vy_std_stripe));
        Bin.vy_prof_af(frame,bin) = nanmean(vy_stripe(vy_inc));
        [mag_inc] = find(mag_stripe >= (mag_mean_stripe-2*mag_std_stripe) & mag_stripe <= (mag_mean_stripe+2*mag_std_stripe));
        Bin.mag_prof_af(frame,bin) = nanmean(mag_stripe(mag_inc));
        
        % Document the number of velocity vectors that were removed in the above steps
        Bin.Remove.vx_stripe(frame,bin) = numel(vx_stripe)-numel(vx_inc);
        Bin.Remove.vy_stripe(frame,bin) = numel(vy_stripe)-numel(vy_inc);
        Bin.Remove.mag_stripe(frame,bin) = numel(mag_stripe)-numel(mag_inc);
            
    end % for m = 1:No_bins
         
end % for Frame

% This final mean will give you the mean velocities of the stripe over time
Bin.vxp_immean = nanmean(Bin.vxp_af);
Bin.vxp_norm_immean = Bin.vxp_immean/sqrt(nanmean(Bin.vxp_af.^2));
Bin.vyp_immean = nanmean(Bin.vyp_af);
Bin.vyp_norm_immean = Bin.vyp_immean/sqrt(nanmean(Bin.vyp_af.^2));
Bin.magp_immean = nanmean(Bin.magp_af);
Bin.magp_chc_immean = nanmean(Bin.magp_chc_af);

Bin.vxa_immean = nanmean(Bin.vxa_af);
Bin.vxa_norm_immean = Bin.vxa_immean/sqrt(nanmean(Bin.vxa_af.^2));
Bin.vya_immean = nanmean(Bin.vya_af);
Bin.vya_norm_immean = Bin.vya_immean/sqrt(nanmean(Bin.vya_af.^2));
Bin.maga_immean = nanmean(Bin.maga_af);

Bin.vxm_immean = nanmean(Bin.vxm_af);
Bin.vxm_norm_immean = Bin.vxm_immean/sqrt(nanmean(Bin.vxm_af.^2));
Bin.vym_immean = nanmean(Bin.vym_af);
Bin.vym_norm_immean = Bin.vym_immean/sqrt(nanmean(Bin.vym_af.^2));
Bin.magm_immean = nanmean(Bin.magm_af);

Bin.vxp_ini_immean = nanmean(Bin.vxp_af(1:10));
Bin.vyp_ini_immean = nanmean(Bin.vyp_af(1:10));
Bin.vxp_end_immean = nanmean(Bin.vxp_af(end-9:end));
Bin.vyp_end_immean = nanmean(Bin.vyp_af(end-9:end));

Bin.vxa_ini_immean = nanmean(Bin.vxa_af(1:10));
Bin.vya_ini_immean = nanmean(Bin.vya_af(1:10));
Bin.vxa_end_immean = nanmean(Bin.vxa_af(end-9:end));
Bin.vya_end_immean = nanmean(Bin.vya_af(end-9:end));

Bin.magp_ini_immean = nanmean(Bin.magp_af(1:10));
Bin.magp_end_immean = nanmean(Bin.magp_af(end-9:end));
Bin.maga_ini_immean = nanmean(Bin.maga_af(1:10));
Bin.maga_end_immean = nanmean(Bin.maga_af(end-9:end));

%% Define x-axis for the graphs
x_axis = 1:No_bins;
x_axis = (x_axis - min(x_axis))/(max(x_axis) - min(x_axis));

clf;plot(x_axis,nanmean(Bin.vx_prof_af),'Color','k','LineWidth',2,'Marker','o','MarkerSize',6,'MarkerEdgeColor',[.2 .2 .2],'MarkerFaceColor',[.7 .7 .7]);
set(gca,'XLim',[-0.2 1.2],'XTick',0:0.2:1,'YDir','reverse')
saveas(gcf,strcat(Folder_path,'/Analysis_v1/VelProf_vx'),'pdf');

clf;plot(x_axis,nanmean(Bin.vy_prof_af),'Color','k','LineWidth',2,'Marker','o','MarkerSize',6,'MarkerEdgeColor',[.2 .2 .2],'MarkerFaceColor',[.7 .7 .7]);
set(gca,'XLim',[-0.2 1.2],'XTick',0:0.2:1)
saveas(gcf,strcat(Folder_path,'/Analysis_v1/Velprof_vy'),'pdf');

close all

save(strcat(Folder_path,'/Analysis_v1/Bin_pivlab.mat'),'Bin')
% disp('------Binning over------')


load(strcat(Folder_path,'/Analysis_v1/Bin_pivlab.mat'))
load(strcat(Folder_path,'/Analysis_v1/Rot_pivlab.mat'))

I = imread(strcat(Folder_path,'/Parameters/A0.tif'));
I = I.*0.6; % This makes the image a little dimmer. I = I.*0.7 for Fig.S3A and I = I.*0.8 for Fig.S3B

[Irot,~,~] = imtransform(I,Rot.Rotation_transformation);
half_stripewidth = 30;

clf;imshow(imcomplement(Irot))
hold on
line(Bin.No_binseq,repmat(Rot.h_top+Rot.ecenter_y-half_stripewidth,1,19),'Color',[247/255 148/255 30/255])
line(Bin.No_binseq,repmat(Rot.h_top+Rot.ecenter_y+half_stripewidth,1,19),'Color',[247/255 148/255 30/255])
hold on
for i = 1:numel(Bin.No_binseq)
    line([Bin.No_binseq(i) Bin.No_binseq(i)],[Rot.h_top+Rot.ecenter_y-half_stripewidth Rot.h_top+Rot.ecenter_y+half_stripewidth],'Color',[247/255 148/255 30/255])
    hold on
end
axis off
set(gcf,'Color','w')
daspect([1 1 1])
saveas(gcf,strcat(Folder_path,'/Analysis_v1/Image_with_bins'),'tif')