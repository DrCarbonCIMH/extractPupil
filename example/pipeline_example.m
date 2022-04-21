function  pipeline_example

% example of a possible pipeline; make sure that the path is set in Matlab
% uses the short video "pupil_example.wmv" and processes it
% writes out diameter, F and resulting video in subfolder "savedir"

% Algorithm by Markus Sack, Robert Becker (IMAGING ZI) ...
% Alternative Algorithm for pupil diameter extraction ... runs maybe more
% stable than the code provided by Matthias Gamer ...

%% PREP ...
% define directories ...
maindir = pwd; % directory where the videos are stored
savedir = [pwd filesep 'savedir']; % directory for saving
% a. get all video files saved in maindir ...
Videolist = getAllFiles(maindir,'*.wmv',1);

%% LOOP OVER VIDEOS STARTS HERE ...
for vid = 1:numel(Videolist)
    %% get video ...
    [Vdir,Vidname,~] = fileparts(Videolist{vid});
    
    %% load the video
    % that is a bit annoying.. even Windows thinks it is 30fps
    v = VideoReader(Videolist{vid});
    
    %% find the pup for the first time to get the center
    video = readFrame(v);
    img = squeeze(video(:,:,1));
    img=imgaussfilt(img,0.5);
    img=imadjust(img);
    s = ms_findPup_rb(img);
    ms_plotPup_rb(img,s,v.CurrentTime,1);
    cent=s.cc;
    % reset to start of video
    v.CurrentTime=0;
    %% START ALGORITHM ...
    tic
    % define a box which is placed around the found center of pupil (could
    % be also s.rch *2 or so; depending on if eye is open or not)
    %     box = 79;
    box = round(s.rch *2);
    % RunTill = end of frame loop ...
    RunTill = v.Duration; % the whole video duration
    ix=1; % indexing
    while v.CurrentTime < RunTill
        % get the next frame
        video = rgb2gray(readFrame(v));
        
        fprintf('Current Time %.2f: ', v.CurrentTime)
        % "box" the image based on the found center
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
            s = ms_findPup_rb(img);
        end
        
        cent=cent+(s.cc-box-1); % update the center
        s.cent=cent; % put it in s for the next iteration
        %% "plot" and save the frame F -> for later video
        F(ix) = ms_plotPup_rb(img,s,v.CurrentTime,1); % plotting to see in "real time"
        %         F(ix)=ms_plotPup_rb(img,s,v.CurrentTime,0); % plotting without actually seeing anything
        radius(ix)=s.rch;
        % step ahead
        ix=ix+1;
    end % end of while v.CurrentTime<=RunTill
    toc
    %% calculate pupil diameter ...
    diameter = radius*2;
    %% save diameter, F, and make a video...
    mkdir(savedir);
    % save the diameter
    save([savedir filesep Vidname '_diameter.mat'],'diameter');
    % save F; maybe F is too big for saving in one mat file ...
    save([savedir filesep Vidname '_Fplay.mat'],'F');
    % make a video
    vidObj = VideoWriter([savedir filesep 'pupil_example.avi']);
    open(vidObj);
    for ix=1:length(F)
        figure(99); imshow(F(ix).cdata);
        currFrame = getframe(gcf);
        writeVideo(vidObj,currFrame);
    end
    
    clear F diameter diameter_lp
end % end video loop ...
