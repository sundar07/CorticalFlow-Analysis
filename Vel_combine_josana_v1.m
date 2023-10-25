function Vel_combine_josana_v1(Gene,No_movies,No_bins)

% function to get a whole lot of statistics from the velocity data
% version 4 - all bin data has been removed as we decided to do analysis only for stripes, and
% have added a few lines for storing individual frame data

%% Pre-allocate the variables

Arby_no_frames = 100;

Vel.post.vx_immean = nan(No_movies,1);
Vel.post.vx_norm_immean = nan(No_movies,1);
Vel.post.vy_immean = nan(No_movies,1);
Vel.post.vy_norm_immean = nan(No_movies,1);
Vel.post.mag_immean = nan(No_movies,1);
Vel.post.mag_chc_immean = nan(No_movies,1);

Vel.ant.vx_immean = nan(No_movies,1);
Vel.ant.vx_norm_immean = nan(No_movies,1);
Vel.ant.vy_immean = nan(No_movies,1);
Vel.ant.vy_norm_immean = nan(No_movies,1);
Vel.ant.mag_immean = nan(No_movies,1);

Vel.mid.vx_immean = nan(No_movies,1);
Vel.mid.vx_norm_immean = nan(No_movies,1);
Vel.mid.vy_immean = nan(No_movies,1);
Vel.mid.vy_norm_immean = nan(No_movies,1);
Vel.mid.mag_immean = nan(No_movies,1);

Vel.post.vx_ini = nan(No_movies,1);
Vel.post.vy_ini = nan(No_movies,1);
Vel.post.vx_end = nan(No_movies,1);
Vel.post.vy_end = nan(No_movies,1);
Vel.ant.vx_ini = nan(No_movies,1);
Vel.ant.vy_ini = nan(No_movies,1);
Vel.ant.vx_end = nan(No_movies,1);
Vel.ant.vy_end = nan(No_movies,1);

Vel.post.vx_af = nan(Arby_no_frames,No_movies);
Vel.post.vy_af = nan(Arby_no_frames,No_movies);
Vel.post.mag_af = nan(Arby_no_frames,No_movies);
Vel.post.mag_chc_af = nan(Arby_no_frames,No_movies);

Vel.ant.vx_af = nan(Arby_no_frames,No_movies);
Vel.ant.vy_af = nan(Arby_no_frames,No_movies);
Vel.ant.mag_af = nan(Arby_no_frames,No_movies);

Vel.mid.vx_af = nan(Arby_no_frames,No_movies);
Vel.mid.vy_af = nan(Arby_no_frames,No_movies);
Vel.mid.mag_af = nan(Arby_no_frames,No_movies);

Vel.Prof.vx_stripe = nan(No_movies,No_bins);
Vel.Prof.vy_stripe = nan(No_movies,No_bins);
Vel.Prof.Mag_stripe = nan(No_movies,No_bins);

Vel.Prof.vx_stripe_af = nan(Arby_no_frames*No_movies,No_bins);
Vel.Prof.vy_stripe_af = nan(Arby_no_frames*No_movies,No_bins);
Vel.Prof.mag_stripe_af = nan(Arby_no_frames*No_movies,No_bins);

%%
RNAi_path = cd;
for movie = 1:No_movies
    
    load(strcat(RNAi_path,'/',Gene,'_',num2str(movie),'/Analysis_v1/Bin_pivlab.mat'))
    
    No_frames = length(Bin.vxp_af);
    
    % Get velocity data from individual frames from each movie in the posterior, middle and anterior
    Vel.post.vx_af(1:No_frames,movie) = Bin.vxp_af;
    Vel.post.vy_af(1:No_frames,movie) = Bin.vyp_af;
    Vel.post.mag_af(1:No_frames,movie) = Bin.magp_af;
    Vel.post.mag_chc_af(1:No_frames,movie) = Bin.magp_chc_af;

    Vel.ant.vx_af(1:No_frames,movie) = Bin.vxa_af;
    Vel.ant.vy_af(1:No_frames,movie) = Bin.vya_af;
    Vel.ant.mag_af(1:No_frames,movie) = Bin.maga_af;

    Vel.mid.vx_af(1:No_frames,movie) = Bin.vxm_af;
    Vel.mid.vy_af(1:No_frames,movie) = Bin.vym_af;
    Vel.mid.mag_af(1:No_frames,movie) = Bin.magm_af;

    % Get mean velocity data from each movie in the posterior and anterior
    Vel.post.vx_immean(movie,1) = Bin.vxp_immean;
    Vel.post.vx_norm_immean(movie,1) = Bin.vxp_norm_immean;
    Vel.post.vy_immean(movie,1) = Bin.vyp_immean;
    Vel.post.vy_norm_immean(movie,1) = Bin.vyp_norm_immean;
    Vel.post.mag_immean(movie,1) = Bin.magp_immean;
    Vel.post.mag_chc_immean(movie,1) = Bin.magp_chc_immean;

    Vel.ant.vx_immean(movie,1) = Bin.vxa_immean;
    Vel.ant.vx_norm_immean(movie,1) = Bin.vxa_norm_immean;
    Vel.ant.vy_immean(movie,1) = Bin.vya_immean;
    Vel.ant.vy_norm_immean(movie,1) = Bin.vya_norm_immean;
    Vel.ant.mag_immean(movie,1) = Bin.maga_immean;

    Vel.mid.vx_immean(movie,1) = Bin.vxm_immean;
    Vel.mid.vx_norm_immean(movie,1) = Bin.vxm_norm_immean;
    Vel.mid.vy_immean(movie,1) = Bin.vym_immean;
    Vel.mid.vy_norm_immean(movie,1) = Bin.vym_norm_immean;
    Vel.mid.mag_immean(movie,1) = Bin.magm_immean;

    Vel.post.vx_ini_immean(movie,1) = Bin.vxp_ini_immean;
    Vel.post.vx_end_immean(movie,1) = Bin.vxp_end_immean;
    Vel.post.vy_ini_immean(movie,1) = Bin.vyp_ini_immean;
    Vel.post.vy_end_immean(movie,1) = Bin.vyp_end_immean;
    Vel.post.mag_ini_immean(movie,1) = Bin.magp_ini_immean;
    Vel.post.mag_end_immean(movie,1) = Bin.magp_end_immean;

    Vel.ant.vx_ini_immean(movie,1) = Bin.vxa_ini_immean;
    Vel.ant.vx_end_immean(movie,1) = Bin.vxa_end_immean;
    Vel.ant.vy_ini_immean(movie,1) = Bin.vya_ini_immean;
    Vel.ant.vy_end_immean(movie,1) = Bin.vya_end_immean;
    Vel.ant.mag_ini_immean(movie,1) = Bin.maga_ini_immean;
    Vel.ant.mag_end_immean(movie,1) = Bin.maga_end_immean;

    % Get mean vel profile data from each movie
    Vel.Prof.vx_prof_immean(movie,:) = nanmean(Bin.vx_prof_af);
    Vel.Prof.vy_prof_immean(movie,:) = nanmean(Bin.vy_prof_af);
    Vel.Prof.mag_prof_immean(movie,:) = nanmean(Bin.mag_prof_af);

    % Get median vel profile data from each movie
    Vel.Prof.vx_prof_immedian(movie,:) = nanmedian(Bin.vx_prof_af);
    Vel.Prof.vy_prof_immedian(movie,:) = nanmedian(Bin.vy_prof_af);
    Vel.Prof.mag_prof_immedian(movie,:) = nanmedian(Bin.mag_prof_af);

    % Get velocity profiles from each frame and all movies
    file_frames = (movie*Arby_no_frames)-(Arby_no_frames-1):(movie*Arby_no_frames)-(Arby_no_frames-1)+(No_frames-1);
    Vel.Prof.vx_prof_af(file_frames,:) = Bin.vx_prof_af;
    Vel.Prof.vy_prof_af(file_frames,:) = Bin.vy_prof_af;
    Vel.Prof.mag_prof_af(file_frames,:) = Bin.mag_prof_af;

    % Get std and sem from each movie for velocity profiles
    Vel.Prof.vx_prof_imstd(movie,:) = nanstd(Bin.vx_prof_af);
    Vel.Prof.vy_prof_imstd(movie,:) = nanstd(Bin.vy_prof_af);
    Vel.Prof.mag_prof_imstd(movie,:) = nanstd(Bin.mag_prof_af);

    Vel.Prof.vx_prof_imsem(movie,:) = Vel.Prof.vx_prof_imstd(movie,:)./sqrt(Find_num_valid(Bin.vx_prof_af));
    Vel.Prof.vy_prof_imsem(movie,:) = Vel.Prof.vy_prof_imstd(movie,:)./sqrt(Find_num_valid(Bin.vy_prof_af));
    Vel.Prof.mag_prof_imsem(movie,:) = Vel.Prof.mag_prof_imstd(movie,:)./sqrt(Find_num_valid(Bin.mag_prof_af));

    % Come out of the folder and go to the next one
    clear Bin
end

%%
% Mean velocities computed from individual movie means for posterior
Vel.post.vx_ammean = nanmean(Vel.post.vx_immean);
Vel.post.vx_amstd = nanstd(Vel.post.vx_immean);
Vel.post.vx_amsem = Vel.post.vx_amstd/sqrt(No_movies);

Vel.post.vx_norm_ammean = nanmean(Vel.post.vx_norm_immean);
Vel.post.vx_norm_amstd = nanstd(Vel.post.vx_norm_immean);
Vel.post.vx_norm_amsem = Vel.post.vx_norm_amstd/sqrt(No_movies);

Vel.post.vy_ammean = nanmean(Vel.post.vy_immean);
Vel.post.vy_amstd = nanstd(Vel.post.vy_immean);
Vel.post.vy_amsem = Vel.post.vy_amstd/sqrt(No_movies);

Vel.post.vy_norm_ammean = nanmean(Vel.post.vy_norm_immean);
Vel.post.vy_norm_amstd = nanstd(Vel.post.vy_norm_immean);
Vel.post.vy_norm_amsem = Vel.post.vy_norm_amstd/sqrt(No_movies);

Vel.post.mag_ammean = nanmean(Vel.post.mag_immean);
Vel.post.mag_amstd = nanstd(Vel.post.mag_immean);
Vel.post.mag_amsem = Vel.post.mag_amstd/sqrt(No_movies);

Vel.post.mag_chc_ammean = nanmean(Vel.post.mag_chc_immean);
Vel.post.mag_chc_amstd = nanstd(Vel.post.mag_chc_immean);
Vel.post.mag_chc_amsem = Vel.post.mag_chc_amstd/sqrt(No_movies);

% Mean velocities computed from individual movie means for posterior_ini
Vel.post.vx_ini_ammean = nanmean(Vel.post.vx_ini_immean);
Vel.post.vx_ini_amstd = nanstd(Vel.post.vx_ini_immean);
Vel.post.vx_ini_amsem = Vel.post.vx_ini_amstd/sqrt(No_movies);

Vel.post.vy_ini_ammean = nanmean(Vel.post.vy_ini_immean);
Vel.post.vy_ini_amstd = nanstd(Vel.post.vy_ini_immean);
Vel.post.vy_ini_amsem = Vel.post.vy_ini_amstd/sqrt(No_movies);

Vel.post.mag_ini_ammean = nanmean(Vel.post.mag_ini_immean);
Vel.post.mag_ini_amstd = nanstd(Vel.post.mag_ini_immean);
Vel.post.mag_ini_amsem = Vel.post.mag_ini_amstd/sqrt(No_movies);

% Mean velocities computed from individual movie means for posterior_end
Vel.post.vx_end_ammean = nanmean(Vel.post.vx_end_immean);
Vel.post.vx_end_amstd = nanstd(Vel.post.vx_end_immean);
Vel.post.vx_end_amsem = Vel.post.vx_end_amstd/sqrt(No_movies);

Vel.post.vy_end_ammean = nanmean(Vel.post.vy_end_immean);
Vel.post.vy_end_amstd = nanstd(Vel.post.vy_end_immean);
Vel.post.vy_end_amsem = Vel.post.vy_end_amstd/sqrt(No_movies);

Vel.post.mag_end_ammean = nanmean(Vel.post.mag_end_immean);
Vel.post.mag_end_amstd = nanstd(Vel.post.mag_end_immean);
Vel.post.mag_end_amsem = Vel.post.mag_end_amstd/sqrt(No_movies);

% Mean velocities computed from individual movie means for anterior_ini
Vel.ant.vx_ini_ammean = nanmean(Vel.ant.vx_ini_immean);
Vel.ant.vx_ini_amstd = nanstd(Vel.ant.vx_ini_immean);
Vel.ant.vx_ini_amsem = Vel.ant.vx_ini_amstd/sqrt(No_movies);

Vel.ant.vy_ini_ammean = nanmean(Vel.ant.vy_ini_immean);
Vel.ant.vy_ini_amstd = nanstd(Vel.ant.vy_ini_immean);
Vel.ant.vy_ini_amsem = Vel.ant.vy_ini_amstd/sqrt(No_movies);

Vel.ant.mag_ini_ammean = nanmean(Vel.ant.mag_ini_immean);
Vel.ant.mag_ini_amstd = nanstd(Vel.ant.mag_ini_immean);
Vel.ant.mag_ini_amsem = Vel.ant.mag_ini_amstd/sqrt(No_movies);

% Mean velocities computed from individual movie means for anterior_end
Vel.ant.vx_end_ammean = nanmean(Vel.ant.vx_end_immean);
Vel.ant.vx_end_amstd = nanstd(Vel.ant.vx_end_immean);
Vel.ant.vx_end_amsem = Vel.ant.vx_end_amstd/sqrt(No_movies);

Vel.ant.vy_end_ammean = nanmean(Vel.ant.vy_end_immean);
Vel.ant.vy_end_amstd = nanstd(Vel.ant.vy_end_immean);
Vel.ant.vy_end_amsem = Vel.ant.vy_end_amstd/sqrt(No_movies);

Vel.ant.mag_end_ammean = nanmean(Vel.ant.mag_end_immean);
Vel.ant.mag_end_amstd = nanstd(Vel.ant.mag_end_immean);
Vel.ant.mag_end_amsem = Vel.ant.mag_end_amstd/sqrt(No_movies);

% Mean velocities computed from individual movie means for anterior
Vel.ant.vx_ammean = nanmean(Vel.ant.vx_immean);
Vel.ant.vx_amstd = nanstd(Vel.ant.vx_immean);
Vel.ant.vx_amsem = Vel.ant.vx_amstd/sqrt(No_movies);

Vel.ant.vx_norm_ammean = nanmean(Vel.ant.vx_norm_immean);
Vel.ant.vx_norm_amstd = nanstd(Vel.ant.vx_norm_immean);
Vel.ant.vx_norm_amsem = Vel.ant.vx_norm_amstd/sqrt(No_movies);

Vel.ant.vy_ammean = nanmean(Vel.ant.vy_immean);
Vel.ant.vy_amstd = nanstd(Vel.ant.vy_immean);
Vel.ant.vy_amsem = Vel.ant.vy_amstd/sqrt(No_movies);

Vel.ant.vy_norm_ammean = nanmean(Vel.ant.vy_norm_immean);
Vel.ant.vy_norm_amstd = nanstd(Vel.ant.vy_norm_immean);
Vel.ant.vy_norm_amsem = Vel.ant.vy_norm_amstd/sqrt(No_movies);

Vel.ant.mag_ammean = nanmean(Vel.ant.mag_immean);
Vel.ant.mag_amstd = nanstd(Vel.ant.mag_immean);
Vel.ant.mag_amsem = Vel.ant.mag_amstd/sqrt(No_movies);

% Mean velocities computed from individual movie means for middle
Vel.mid.vx_ammean = nanmean(Vel.mid.vx_immean);
Vel.mid.vx_amstd = nanstd(Vel.mid.vx_immean);
Vel.mid.vx_amsem = Vel.mid.vx_amstd/sqrt(No_movies);

Vel.mid.vx_norm_ammean = nanmean(Vel.mid.vx_norm_immean);
Vel.mid.vx_norm_amstd = nanstd(Vel.mid.vx_norm_immean);
Vel.mid.vx_norm_amsem = Vel.mid.vx_norm_amstd/sqrt(No_movies);

Vel.mid.vy_ammean = nanmean(Vel.mid.vy_immean);
Vel.mid.vy_amstd = nanstd(Vel.mid.vy_immean);
Vel.mid.vy_amsem = Vel.mid.vy_amstd/sqrt(No_movies);

Vel.mid.vy_norm_ammean = nanmean(Vel.mid.vy_norm_immean);
Vel.mid.vy_norm_amstd = nanstd(Vel.mid.vy_norm_immean);
Vel.mid.vy_norm_amsem = Vel.mid.vy_norm_amstd/sqrt(No_movies);

Vel.mid.mag_ammean = nanmean(Vel.mid.mag_immean);
Vel.mid.mag_amstd = nanstd(Vel.mid.mag_immean);
Vel.mid.mag_amsem = Vel.mid.mag_amstd/sqrt(No_movies);

%% statistics for averages in stripes

Vel.Prof.vx_prof_ammean = nanmean(Vel.Prof.vx_prof_immean);
Vel.Prof.vx_prof_amstd = nanstd(Vel.Prof.vx_prof_immean);
for i = 1:No_bins
    Vel.Prof.vx_prof_amsem(1,i) = Vel.Prof.vx_prof_amstd(1,i)/sqrt(numel(find(~isnan(Vel.Prof.vx_prof_immean(:,i)))));
end

Vel.Prof.vy_prof_ammean = nanmean(Vel.Prof.vy_prof_immean);
Vel.Prof.vy_prof_amstd = nanstd(Vel.Prof.vy_prof_immean);
for i = 1:No_bins
    Vel.Prof.vy_prof_amsem(1,i) = Vel.Prof.vy_prof_amstd(1,i)/sqrt(numel(find(~isnan(Vel.Prof.vy_prof_immean(:,i)))));
end

Vel.Prof.mag_prof_ammean = nanmean(Vel.Prof.mag_prof_immean);
Vel.Prof.mag_prof_amstd = nanstd(Vel.Prof.mag_prof_immean);
for i = 1:No_bins
    Vel.Prof.mag_prof_amsem(1,i) = Vel.Prof.mag_prof_amstd(1,i)/sqrt(numel(find(~isnan(Vel.Prof.mag_prof_immean(:,i)))));
end


save(strcat(RNAi_path,'/Analysis_v1/Vel_pivlab.mat'),'Vel')

%% Graph time :)

x = normalize(1:No_bins);

% Velocity profiles (ammean)
clf;errorbar(x,Vel.Prof.vx_prof_ammean,Vel.Prof.vx_prof_amsem,'Color','k','LineWidth',2,'Marker','o','MarkerSize',6,'MarkerEdgeColor',[.2 .2 .2],'MarkerFaceColor',[.7 .7 .7])
set(gca,'XLim',[-0.2 1.2],'XTick',0:0.2:1,'YDir','reverse')
saveas(gcf,strcat(RNAi_path,'/Analysis_v1/VelProf_vx_ammean'),'pdf')

clf;errorbar(x,Vel.Prof.vy_prof_ammean,Vel.Prof.vy_prof_amsem,'LineWidth',2,'Marker','o','MarkerSize',6,'MarkerEdgeColor',[.2 .2 .2],'MarkerFaceColor',[.7 .7 .7])
set(gca,'XLim',[-0.2 1.2],'XTick',0:0.2:1)
saveas(gcf,strcat(RNAi_path,'/Analysis_v1/VelProf_vy_ammean'),'pdf')

clf;errorbar(x,Vel.Prof.mag_prof_ammean,Vel.Prof.mag_prof_amsem,'LineWidth',2,'Marker','o','MarkerSize',6,'MarkerEdgeColor',[.2 .2 .2],'MarkerFaceColor',[.7 .7 .7])
set(gca,'XLim',[-0.2 1.2],'XTick',0:0.2:1)
saveas(gcf,strcat(RNAi_path,'/Analysis_v1/VelProf_mag_ammean'),'pdf')

close all

% Make histograms of anterior and posterior AP and orthogonal velocities
% data1 = Vel.ant.vx_af(:);
% data2 = Vel.post.vx_af(:);
% 
% Mi = floor(min(min(data1),min(data2)));
% Ma = floor(max(max(data1),max(data2)));
% 
% x = Mi:(Ma-Mi)/12:Ma;
% 
% norm1 = normpdf(x,nanmean(data1),nanstd(data1));
% nf1 = sum(norm1);
% norm1 = normpdf(Mi:0.0001:Ma,nanmean(data1),nanstd(data1));
% norm1 =norm1./nf1;
% 
% n1 = hist(data1,x);
% n1 = n1./sum(n1);
% 
% norm2 = normpdf(x,nanmean(data2),nanstd(data2));
% nf2 = sum(norm2);
% norm2 = normpdf(Mi:0.0001:Ma,nanmean(data2),nanstd(data2));
% norm2 = norm2./nf2;
% 
% n2 = hist(data2,x);
% n2 = n2./sum(n2);
% 
% clf;
% h1 = bar(x,n1,1,'edgecolor',[0.1 0.1 0.1],'facecolor',[102/255 195/255 162/255]);
% hold on
% plot(Mi:0.0001:Ma,norm1,'Color',[102/255 195/255 162/255])
% hold on
% h2 = bar(x,n2,1,'edgecolor',[0.1 0.1 0.1],'facecolor',[241/255 91/255 106/255]);
% hold on
% plot(Mi:0.0001:Ma,norm2,'Color',[241/255 91/255 106/255])
% 
% alpha(get(h1,'children'),.25)
% alpha(get(h2,'children'),.25)
% set(gcf,'Color','w')
% set(gca,'box','off')
% set(gcf,'Position',[88 264 375 375])
% export_fig(strcat(RNAi_path,'/Analysis_v1/vx_post_ant_hist.pdf'),'-q101')
% 
% clear data1 data2
% 
% data1 = Vel.ant.vy_af(:);
% data2 = Vel.post.vy_af(:);
% 
% Mi = floor(min(min(data1),min(data2)));
% Ma = floor(max(max(data1),max(data2)));
% 
% x = Mi:(Ma-Mi)/12:Ma;
% 
% norm1 = normpdf(x,nanmean(data1),nanstd(data1));
% nf1 = sum(norm1);
% norm1 = normpdf(Mi:0.0001:Ma,nanmean(data1),nanstd(data1));
% norm1 =norm1./nf1;
% 
% n1 = hist(data1,x);
% n1 = n1./sum(n1);
% 
% norm2 = normpdf(x,nanmean(data2),nanstd(data2));
% nf2 = sum(norm2);
% norm2 = normpdf(Mi:0.0001:Ma,nanmean(data2),nanstd(data2));
% norm2 = norm2./nf2;
% 
% n2 = hist(data2,x);
% n2 = n2./sum(n2);
% 
% clf;
% h1 = bar(x,n1,1,'edgecolor',[0.1 0.1 0.1],'facecolor',[102/255 195/255 162/255]);
% hold on
% plot(Mi:0.0001:Ma,norm1,'Color',[102/255 195/255 162/255])
% hold on
% h2 = bar(x,n2,1,'edgecolor',[0.1 0.1 0.1],'facecolor',[241/255 91/255 106/255]);
% hold on
% plot(Mi:0.0001:Ma,norm2,'Color',[241/255 91/255 106/255])
% 
% alpha(get(h1,'children'),.25)
% alpha(get(h2,'children'),.25)
% set(gcf,'Color','w')
% set(gca,'box','off')
% set(gcf,'Position',[88 264 375 375])
% axis square
% 
% export_fig(strcat(RNAi_path,'/Analysis_v1/vy_post_ant_hist.pdf'),'-q101')
% 
% clear data
% 
% Make boxplots of anterior and posterior velocities at different time points
% data(:,1) = Vel.ant.vx_immean;
% data(:,2) = Vel.ant.vx_ini_immean;
% data(:,3) = Vel.ant.vx_end_immean;
% data(:,4) = Vel.post.vx_immean;
% data(:,5) = Vel.post.vx_ini_immean;
% data(:,6) = Vel.post.vx_end_immean;
% clf;notBoxPlot(data,[1 1.75 2.5 3.25 4 4.75])
% set(gcf,'Position',[88 264 375 375])
% axis square
% set(gcf,'Color','w')
% set(gca,'XLim',[0.5 5.25],'XTick',[])
% export_fig(strcat(RNAi_path,'/Analysis_v1/vx_post_ant_ini_end_full.pdf'),'-q101')
% 
% clear data
% 
% data(:,1) = Vel.ant.vy_immean;
% data(:,2) = Vel.ant.vy_ini_immean;
% data(:,3) = Vel.ant.vy_end_immean;
% data(:,4) = Vel.post.vy_immean;
% data(:,5) = Vel.post.vy_ini_immean;
% data(:,6) = Vel.post.vy_end_immean;
% clf;notBoxPlot(data,[1 1.75 2.5 3.25 4 4.75])
% set(gcf,'Position',[88 264 375 375])
% axis square
% set(gcf,'Color','w')
% set(gca,'XLim',[0.5 5.25],'XTick',[])
% export_fig(strcat(RNAi_path,'/Analysis_v1/vy_post_ant_ini_end_full.pdf'),'-q101')
% 
% clear data

% disp('------Vel Analysis_v1 done------')
