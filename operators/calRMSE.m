function [ y ] = calRMSE( a,b )
% function y = calRMSE(a,b)
%
% calculate the RMSE of a and b


if ~isequal(size(a),size(b))
    y = -1;
else
    if length(size(a)) == 2 || length(size(a)) == 1
        diff = abs(a-b).^2;
        y = sqrt(mean(diff(:)));
    end
    if length(size(a)) == 3
        for im = 1:size(a,3)
            diff = abs(squeeze(a(:,:,im))-squeeze(b(:,:,im))).^2;
            y(im) = sqrt(mean(diff(:)));
        end
    end
end
