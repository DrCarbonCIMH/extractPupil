function [img] = ms_preprocIMG(img,s)

if(nargin>1)
    vig = ms_2Dgaussian(floor(size(img,1)/2), [1 1], s);
    vig = -vig./max(vig(:))+2;
    img=cast(img,'single').*vig;
    img = cast(img,'uint8');
end

img=medfilt2(img);
img=imadjust(img);