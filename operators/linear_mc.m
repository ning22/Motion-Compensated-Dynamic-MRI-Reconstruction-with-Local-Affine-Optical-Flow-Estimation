function I1 = linear_mc(I0,hor_ind, ver_ind, ROW, COL)

% global ROW COL
% [ROW COL] = size(I0);

% Interpolated locations in vectorized form
ver_floor = min(floor(ver_ind(:)),ROW-1);
ver_ceil = ver_floor+1;
hor_floor = min(floor(hor_ind(:)),COL-1);
hor_ceil = hor_floor+1;

% Consider a box of pixels surrounding the interpolating pixel
int_ind1 = ver_floor + ROW*(hor_floor-1); % top left corner
int_ind2 = ver_ceil + ROW*(hor_floor-1); % bottom left corner
int_ind3 = ver_floor + ROW*(hor_ceil-1); % top right corner
int_ind4 = ver_ceil + ROW*(hor_ceil-1); % bottom right corner

% Column interpolation
% 1D --> x(ai) = x(ai_floor)*(ai_ceil-ai) + x(ai_ceil)*(ai-ai_floor);
I1_left = I0(int_ind1).*(ver_ceil-ver_ind(:)) + I0(int_ind2).*(ver_ind(:)-ver_floor);
I1_right = I0(int_ind3).*(ver_ceil-ver_ind(:)) + I0(int_ind4).*(ver_ind(:)-ver_floor);

% Row interpolation
% 1D --> x(ai) = x(ai_floor)*(ai_ceil-ai) + x(ai_ceil)*(ai-ai_floor);
I1_temp = I1_left.*(hor_ceil-hor_ind(:)) + I1_right.*(hor_ind(:)-hor_floor);
I1 = reshape(I1_temp,ROW,COL);
