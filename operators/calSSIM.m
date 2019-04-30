function [ ssim ] = calSSIM( image1,image2 )
% function ssim = calSSIM(image1,image2)
% 
% calculate the mean structure similarity index (ssim) between two images
% and two image series


if size(image1) ~= size(image2),
    ssim = -1;
else
    image1 = abs(image1);
    image2 = abs(image2);
    if length(size(image1)) ~=3       
        im1 = (image1-min(image1(:)))*255/(max(image1(:))-min(image1(:)));
        im2 = (image2-min(image2(:)))*255/(max(image2(:))-min(image2(:)));
        ssim = ssim_index(im1,im2);
    end
    if length(size(image1)) == 3,
        nim = size(image1,3);
        for im = 1:nim
            im1 = squeeze(image1(:,:,im));
            im1 = (im1-min(im1(:)))*255/(max(im1(:))-min(im1(:)));
            im2 = squeeze(image2(:,:,im));
            im2 = (im2-min(im2(:)))*255/(max(im2(:))-min(im2(:)));
            ssim(im) = ssim_index(im1,im2);
        end
    end
end

