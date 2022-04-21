function  ExtractPupil_LW_4NatCom

% March 2019 by LW

% Algorithm by Markus Sack, Robert Becker (IMAGING ZI) ...
% Alternative Algorithm for pupil diameter extraction ... runs maybe more
% stable than the code provided by Matthias Gamer ...


%%
clear all;

%% PREP ...

% define directories ...
maindir = '\\zi\flstorage\dep_psychiatrie_psychotherapie\group_entwbio\data_ZIFNAS\Laurens W\Analysis_PupilDiameter\ZI_algorithm\videodir'; % all videos in maindir get processed ...
savedir = '\\zi\flstorage\dep_psychiatrie_psychotherapie\group_entwbio\data_ZIFNAS\Laurens W\Analysis_PupilDiameter\ZI_algorithm\results'
maindir = pwd; savedir = [pwd filesep 'savedir'];
% a. get all video files saved in maindir ...
Videolist = getAllFiles(maindir,'rw*.wmv',1);


%% LOOP OVER VIDEOS STARTS HERE ...
for vid = 1:numel(Videolist)
    %% get general info about current video ...
    [Vdir,Vidname,~] = fileparts(Videolist{vid});
    VideoPathFull = Videolist{vid};
    find_ = strfind(Vidname,'_');
    subjcur = Vidname(1:find_(1)-1);
    date = Vidname(find_(1)+1:end);
    
    
    %% load the video
    % that is a bit annoying.. even Windows thinks it is 30fps
    v=VideoReader(Videolist{vid});
    % you can convert the video to an avi file in VLC
    % v=VideoReader('rw2_181012.avi');
    
    
    %% find the pup for the first time to get the center
    % v.CurrentTime=0;
    video = readFrame(v);
    img = squeeze(video(:,:,1));
    img=imgaussfilt(img,0.5);
    img=imadjust(img);
    s = ms_findPup_rb(img);
    ms_plotPup_rb(img,s,v.CurrentTime,1);
    cent=s.cc;
    v.CurrentTime=0;
    
    
    
    %% START ALGORITHM ...
    
    tic
    
    figure(21); close(21); figure(21); ax=gca; clear pup T F
    % define a box which is placed around the found center of pupil (could
    % be also s.rch *2 or so; depending on if eye is open or not)
    box = 79;
    ix=1;
    
    
    % RunTill = end of frame loop ...
    RunTill = v.Duration-0.05; % remove last frame recorded -> avoid index problems ...
    
    % v.CurrentTime = v.Duration - 2;
    
    while v.CurrentTime<=RunTill
        
        video = rgb2gray(readFrame(v));
        
        fprintf('Current Time %.2f: ', v.CurrentTime)
        T(ix) = v.CurrentTime;
        img = squeeze(video(floor(cent(2)-box):floor(cent(2)+box),floor(cent(1)-box):floor(cent(1)+box),1));
        
        % image adjustments
        [img] = ms_preprocIMG(img,s.rch*1.4);
        
        % actual pup finding
        if(exist('s', 'var'))
            if(~isfield('s','fallBack'))
                s = ms_findPup_rb(img,s);
            else
                s = ms_findPup_rb(img);
            end
        else
            s = ms_findPup(img);
        end
        
        cent=cent+(s.cc-box-1); % update the center
        s.cent=cent; % needed?
        
        
        %% live view
        %     pup(ix,1)=s.Area;
        %     pup(ix,2)=s.MajorAxisLength;
        %     pup(ix,3)=s.ConvexArea;
        %     pup(ix,4)=s.EllArea;
        %     pup(ix,5)=pi*(s.MajorAxisLength/2)^2;
        %     plot(ax,T,[pup(:,1)/pup(1,1), pup(:,2)/pup(1,2), pup(:,3)/pup(1,3), pup(:,4)/pup(1,4), pup(:,5)/pup(1,5)]); %legend({'Area', 'MajorAxis'})
        %
        %     pup(ix,1)=s.rch;
        %     pup(ix,2)=2*pi*(s.rch).^2;
        %     pup(ix,3)=s.Area;
        %     plot(ax,T,[pup(:,1)/pup(1,1), pup(:,2)/pup(1,2), pup(:,3)/pup(1,3)]); %legend({'Convex Hull Radius','Convex Hull Area','Area'});
        
        %%
        F(ix)=ms_plotPup_rb(img,s,v.CurrentTime,0); % plotting
        
        
        radius(ix)=s.rch;
        % for not live videos
        %     s.img = img;
        %     pup{ix} =s;
        % step ahead
        ix=ix+1;
        
        % jump in time
        %      v.CurrentTime=(ix-1)*1;
    end % end of while v.CurrentTime<=RunTill
    toc
    
    
    %% calculate pupil diameter ...
    diameter = radius*2;
    
    %% filtering diameter ...
    h = fdesign.lowpass('N,F3dB',12,0.15);
    d1 = design(h,'butter');
    diameter_lp = filtfilt(d1.sosMatrix,d1.ScaleValues,diameter); % low pass filter ...
    
    
    %% save diameter, F ...
    
    subjdir = [savedir filesep Vidname];
    mkdir(subjdir);
    
    save([savedir filesep subjdir filesep Vidname '_diameter.mat'],'diameter','diameter_lp');
    
    % save F; maybe F is too big for saving in one mat file ...
    save([savedir filesep subjdir filesep Vidname '_Fplay.mat'],'F');
    % first minute ...
    F_first = F(1:1:(20*1*60));
    save([savedir filesep subjdir filesep Vidname '_FplayFIRSTMIN.mat'],'F_first');
    
    % middle ... (minute 15)
    F_middle = F((20*15*60):1:(20*16*60));
    save([savedir filesep subjdir filesep Vidname '_FplayMIDDLEMIN.mat'],'F_middle');
    
    % last minute ...
    F_end = F((end - (20*1*60)):1:end);
    save([savedir filesep subjdir filesep Vidname '_FplayLASTMIN.mat'],'F_end');
    
    
    
    clear F diameter diameter_lp
    
    
    
end % end video loop ...

