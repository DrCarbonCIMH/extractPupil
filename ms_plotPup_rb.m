function [F] = ms_plotPup(bw,s,t, show)
%MS_PLOTPUP Summary of this function goes here
%   Detailed explanation goes here
if show
figure(23)
else
    h=figure('Visible', 'off');
end
imshow([bw (s.PCNN)*255],'InitialMagnification','fit')
% imshowpair([bw s.PCNN],video,'montage')
% title(sprintf('Area %0.2f', s.Area))
text(0,10,sprintf('T: %.2f; Area: %.0f; ',t, s.Area),'Color','red','FontSize',14)
text(0,28,sprintf('R: %.2f;', s.rch),'Color','red','FontSize',14)
% t = linspace(0,2*pi,50);

hold on
for k = 1:length(s)
%     a = s(k).MajorAxisLength/2;
%     b = s(k).MinorAxisLength/2;
%     Xc = s(k).Centroid(1);
%     Yc = s(k).Centroid(2);
%     phi = deg2rad(-s(k).Orientation);
%     x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
%     y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
%     plot(x,y,'r','Linewidth',2)
%     plot(s.Centroid(1),s.Centroid(2),'o')
%     MAx = 0; MAy = -a:1:a;
%     [THETA,R] = cart2pol(MAx,MAy); %Convert to polar coordinates
% %     THETA=THETA-phi-pi/2; %Add a_rad to theta
%     THETA=THETA+phi+pi/2; %Add a_rad to theta
%     [xr,yr] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
%     plot(xr+Xc,yr+Yc,'g-'); hold on; %Original

    viscircles(s(k).cc,s(k).rch,'color','r');

end
hold off

% imshow(bw)
% hold on
% 
% phi = linspace(0,2*pi,50);
% cosphi = cos(phi);
% sinphi = sin(phi);
% 
% for k = 1:length(s)
%     xbar = s(k).Centroid(1);
%     ybar = s(k).Centroid(2);
% 
%     a = s(k).MajorAxisLength/2;
%     b = s(k).MinorAxisLength/2;
% 
%     theta = pi*s(k).Orientation/180;
%     R = [ cos(theta)   sin(theta)
%          -sin(theta)   cos(theta)];
% 
%     xy = [a*cosphi; b*sinphi];
%     xy = R*xy;
% 
%     x = xy(1,:) + xbar;
%     y = xy(2,:) + ybar;
% 
%     plot(x,y,'r','LineWidth',2);
% end
% hold off

F = getframe;
if ~show; close(h); end;

end

