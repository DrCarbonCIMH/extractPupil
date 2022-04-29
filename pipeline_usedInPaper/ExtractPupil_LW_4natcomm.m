edit function  ExtractPupil_ZIalgorithm_LW

% Algorithm by Markus Sack, Robert Becker ... 
% Algorithm for pupil diameter extraction 


%% 
clear all; 

%% PREP ... 

% define directories ... 
maindir = '\\zi\flstorage\dep_psychiatrie_psychotherapie\group_entwbio\data_ZIFNAS\Laurens W\Analysis_PupilDiameter\ZI_algorithm\videodir'; % all videos in maindir get processed ... 
savedir = '\\zi\flstorage\dep_psychiatrie_psychotherapie\group_entwbio\data_ZIFNAS\Laurens W\Analysis_PupilDiameter\ZI_algorithm\results'

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
v=VideoReader(Videolist{vid});


%% find the pup for the first time to get the center
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
box = 79;
ix=1;


% RunTill = end of frame loop ... 
RunTill = v.Duration-0.05; % remove last frame recorded -> avoid index problems ... 


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
    ix=ix+1;
    
end
toc


%% calculate pupil diameter ... 
diameter = radius*2; 


 end % end video loop ... 

