function Master_analysis_v1()

% This function does all the things given by the following m-files and this
% function is used to analyze multiple movies stored in different folders

Folder_path = cd;

DirOutput = dir(strcat(Folder_path,'/*_rnai'));
RNAi_Names = {DirOutput.name}';
Num_Knockdowns = numel(RNAi_Names);

half_stripewidth = 30; % pixels
pixelsize = 0.11;
frameinterval = 1;
No_bins = 18;

for RNAi = 1:Num_Knockdowns
    
    Gene = strrep(RNAi_Names{RNAi},'_rnai','');
    
    cd(strcat(Folder_path,'/',RNAi_Names{RNAi}))
    RNAi_path = cd;
    
    Movie_folders = dir(strcat(RNAi_path,'/',Gene,'_*'));
    Movie_foldername = {Movie_folders.name};
    
    for Movies = 1:numel(Movie_foldername)
        name = Movie_foldername{Movies};
        a = strfind(name,'_');
        Num_underscores = length(a);
        b(Movies,1) = isempty(str2num(name(a(Num_underscores)+1:end)));
    end
    
    No_movies = numel(Movie_foldername)-sum(b);
    
    for movie = 1:No_movies

        cd(strcat(RNAi_path,'/',Gene,'_',num2str(movie)))
        movie_path = cd;

        if exist(strcat(movie_path,'/Analysis_v1'),'dir') ~= 7
           mkdir(strcat(movie_path,'/Analysis_v1'))
        end

        % Open the text file named 'parameters.txt' and read in the parameters
        file = fopen(strcat(movie_path,'/Parameters/parameters.txt'));
        degree = double(cell2mat(textscan(file,'%*s %d',1)));
        % frames = double(cell2mat(textscan(file,'%*s %*s %d %d','HeaderLines',1)));
        fseek(file,0,'bof');
        % frameinterval = double(cell2mat(textscan(file,'%*s %d','HeaderLines',2)));
        fclose(file);
        
        specifyframes = 1:100;
        % specifyframes = frames(1):frames(2);
        
        PIVlab_commandline_josana_v1(specifyframes)
        load(strcat(movie_path,'/Analysis_v1/Results_PIVlab/PIVlab_res.mat'))

        for i = 1:length(specifyframes)-1
            x(:,:,i) = PIVlab_res.x{i};
            y(:,:,i) = PIVlab_res.y{i};
            u(:,:,i) = PIVlab_res.u_filt{i};
            v(:,:,i) = PIVlab_res.v_filt{i};
        end

        Rot_vec_josana_v1(x,y,u,v,degree,PIVlab_res.mask)
        load(strcat(movie_path,'/Analysis_v1/Rot_pivlab.mat'))

        Bin_vel_josana_v1(Rot.xrot,Rot.yrot,Rot.vxrot,Rot.vyrot,No_bins,specifyframes,Rot.Ppole,Rot.Apole,Rot.h_top,Rot.h_bottom,pixelsize,frameinterval,half_stripewidth)
        
        clear PIVlab_res Rot_pivlab Bin_pivlab x y u v
        cd ..
        cprintf('comment',['End of movie ' num2str(movie) ' (Total - ' num2str(No_movies) 'movies)\n']);
    end % for movie

    if exist(strcat(RNAi_path,'/Analysis_v1'),'dir') ~= 7
       mkdir(strcat(RNAi_path,'/Analysis_v1'))
    end

    Vel_combine_josana_v1(Gene,No_movies,No_bins)
    
    clear Vel b
    cprintf('comment',['End of ' RNAi_Names{RNAi} ' (' num2str(RNAi) '/' num2str(Num_Knockdowns) ' knockdowns)\n']);
    cd ..
end