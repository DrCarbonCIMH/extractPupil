function [s_old] = ms_findPup(img, Spre)
%MS_FINDPUP Summary of this function goes here
%   Detailed explanation goes here
verb=1; ToSkip=1; 
% S=cast(imcomplement(img(100:end-100,100:end-100)),'single');
S=cast(imcomplement(img),'single');
S=S./max(max(S));
voxdim=size(S); r=1; totVol=prod(voxdim);
N=60; % doesn't look like we need more than 60 iterations
M = ms_2Dgaussian(r, voxdim./min(voxdim), 0.5);

% Y=zeros(size(S)); F = Y; L = F;  T = Y + 1; %U = F;
% A = logical(Y); BW=false([N, size(S,1), size(S,2)]);
% t =1; brainSize=0;
% % for t = 1:N
% % while brainSize<10000
% found=0;
% % s=struct;
% while ~found
%     [F,L,T,Y] = calcPCNN(F,L,T,Y,M,S);
%     A = A | logical(Y);
%     if t>ToSkip % to save time we could skip the first iterations
%         [~, s] = calcBrainMask(A); % this takes the most time ~1-2s for high res images
%         if(~isempty(s))
%             if(exist('s_old','var'))
%                 if(t>ToSkip+1) 
%                     if((s_old.Area/s_old.ConvexArea > s.Area/s.ConvexArea && s_old.Area>1500) || s.Area>10000)
%                         fpfrintf('#%.0f:relation Area/FilledArea %0.3f; Area %.0f FOUND!\n',t, s_old.Area/s_old.ConvexArea, s_old.Area)
%                         break; 
%                     end
%                 end
%             end
% %             fprintf('#%.0f:relation Area/FilledArea %0.3f; Area %.0f\n',t, s.Area/s.ConvexArea, s.Area)
%             s_old=s;
%         end
% %         [BW(t+1,:,:)] = imdilate(imerode(A, strel(6)), strel(6));
% %         vox(t+1)=length(find(BW(t+1,:,:))); %brainSize=vox(t+1);
% %         vox(t+1)=s.Area; brainSize = s.Area;
% %         fprintf('#%.0f: BrainSize: %.2f mm^3 (%.2f of total Volume)\n',t, vox(t+1), vox(t+1)/totVol)
%         if verb; plotSomething(A, A, Y, T,F,L,1,t,S); end
% %         if vox(t+1)>BrSize(2)/1000; break; end
% %         if t>30 && vox(t)>0
% %             if vox(t)<=vox(t-10)*1.01
% %                 vox(end)=BrSize(2)/1000*1.1; warndlg(sprintf('I was trapped and stopped the iteration!\n%s',name)); break; 
% %             end
% %         end
%     end
%     t=t+1;
% %     pause(1)
% end

%% another approach
% another approach could be to simply run the PCNN till the "best
% candidate" area is > 10000 while collecting all bestCandidates and then
% decide which is really the best
Y=zeros(size(S)); F = Y; L = F;  T = Y + 1; %U = F;
A = logical(Y);
area=0; t=1;
M=fspecial('gaussian',13,1);
% while area < 10000
% we could really need a better "break" condition!
% can we roughly estimate the pup size?
re = sum(S(S>0.80))*1.5; fprintf("Esitmate: %.0f ",re)
if(nargin==1) % that means we have no idea
    while sum(A(:))/numel(A) < 0.25
        [F,L,T,Y] = calcPCNN(F,L,T,Y,M,S);
        A = A | logical(Y); %figure; imagesc(A)
        s = bestCandidate(A);
        if(~isempty(s))
            %         area = s.Area;
            s.cc=mean([max(s.ConvexHull);min(s.ConvexHull)]);
            s.dch=((s.ConvexHull(:,1)-s.cc(1)).^2+(s.ConvexHull(:,2)-s.cc(2)).^2).^(1/2);
            s.vdch=var(s.dch);
            s.rch=mean(s.dch);
            
            cands(t)=s;
            t=t+1;
            %         if(length(cands)>1)
            %             if(cands(end-1).Area/cands(end-1).ConvexArea > cands(end).Area/cands(end).ConvexArea && cands(end-1).Area>1500)
            %                 break;
            %             end
            %         end
        end
    end
else
%     while sum(A(:)) < Spre.Area*2 && sum(A(:))/numel(A) < 0.5
    while sum(A(:)) < re
        [F,L,T,Y] = calcPCNN(F,L,T,Y,M,S);
        A = A | logical(Y);
        s = bestCandidate(A);
        if(~isempty(s))
            s.cc=mean([max(s.ConvexHull);min(s.ConvexHull)]);
            s.dch=((s.ConvexHull(:,1)-s.cc(1)).^2+(s.ConvexHull(:,2)-s.cc(2)).^2).^(1/2);
            s.vdch=var(s.dch);
            s.rch=mean(s.dch);
            
            cands(t)=s;
            t=t+1;
        end
    end
end

% get the best of the cands
if(exist('cands','var')); cands(end)=[]; else; fallBack = Spre; fallBack.fallBack=1; fprintf('Could not find anything!\n'); s_old=fallBack; return; end % the last one should be non-sense anyway
if(~isempty(s)); fallBack = s; else; fallBack = Spre; end
fallBack.Centroid = size(A)/2;
fallBack.MajorAxisLength= size(A,1)/2;
fallBack.MinorAxisLength= size(A,1)/2;
fallBack.fallBack=1;
if(isempty(cands)); fprintf('Could not find anything!\n'); s_old=fallBack; return; end
% the following limits could also be dependend on the "result" before
if(nargin==1)
    % if it's too big
    cands([cands.Area]>6000) = [];
else
    % if it's too far away
    for ix=1:length(cands); xr(ix)=cands(ix).cc(1); yr(ix)=cands(ix).cc(2); end
    cands(xr>150 | xr<50 | yr>150 | yr<50) = [];
%     cands([cands.MajorAxisLength]>Spre.MajorAxisLength*1.5) = [];
end
% relation of Area and ConvexArea seems to be a good "measure"
% [~,idx]=sort([cands.Area]./[cands.ConvexArea]);

% eccentricity may be better -RB
% [~,idx]=sort([cands.Eccentricity],'descend');

% use "roundness" of convex hull -RB
[~,idx]=sort([cands.vdch],'descend');

try
    s_old=cands(idx(end));
    fprintf(' Radius %.0f CHOSEN#%.0f of %.0f\n', s_old.rch,idx(end), t)
catch % as a fallback
    s_old=fallBack;
    fprintf('Could not find anything!\n');
end

end

function s = bestCandidate(A)
% A = bwmorph(A,'clean');
% A = bwmorph(A,'spur');
% A = bwmorph(A,'hbreak');

se = strel('cube',3);
A = imdilate(imerode(A,se), se);
A = imfill(A,'holes');
% s = regionprops(A,{...
%     'Centroid',...
%     'MajorAxisLength',...
%     'MinorAxisLength',...
%     'Orientation', 'Area', 'ConvexArea', 'ConvexImage', 'Eccentricity','ConvexHull'});

s = regionprops(A,{'Area','ConvexArea','ConvexHull'});

% quick but a bit dirty
% s = regionprops(A,{...
%     'Centroid',...
%     'MajorAxisLength',...
%     'MinorAxisLength',...
%     'Orientation', 'Area'});
% for ix=1:length(s); s(ix).ConvexArea = s(ix).MajorAxisLength*s(ix).MinorAxisLength*pi/4; end

% s=s([s.FilledArea]==max([s.FilledArea]));
% s([s.MajorAxisLength]<40)=[]; % remove if too small
s([s.Area]<600)=[]; % remove if too small
% if(sum([s.MajorAxisLength]./[s.MinorAxisLength]<1.1)>0); s([s.MajorAxisLength]./[s.MinorAxisLength]>1.1)=[]; end% remove if not somwhat a circle
% best guess seems to be the closest relation of Area/ConvexArea to 1
s=s([s.Area]./[s.ConvexArea]==max([s.Area]./[s.ConvexArea]));
s=s([s.Area] == max([s.Area]));

if(~isempty(s)); s.PCNN=A; end % s.EllArea = s.MajorAxisLength*s.MinorAxisLength*pi/4; end
end

