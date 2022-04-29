function [F] = ms_plotPup_rb(bw,s,t, show)

if show
figure(23)
else
    h=figure('Visible', 'off');
end
imshow([bw (s.PCNN)*255],'InitialMagnification','fit')
text(0,10,sprintf('T: %.2f; Area: %.0f; ',t, s.Area),'Color','red','FontSize',14)
text(0,28,sprintf('R: %.2f;', s.rch),'Color','red','FontSize',14)

hold on
for k = 1:length(s)
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