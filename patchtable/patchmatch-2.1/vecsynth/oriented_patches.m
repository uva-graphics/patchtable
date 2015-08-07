function patches = oriented_patches(I, angle, patch_w, filter_mode)

h = size(I, 1);
w = size(I, 2);

patch_dims = patch_w^2;

eh = h-patch_w+1;
ew = w-patch_w+1;

hw = floor(patch_w/2);

patches = zeros(eh, ew, patch_dims);
angle = angle*pi/180;
cos_angle = cos(angle);
sin_angle = sin(angle);
[du_grid, dv_grid] = meshgrid((-hw):hw, (-hw):hw);
du_grid = du_grid(:);
dv_grid = dv_grid(:);

if ~exist('filter_mode', 'var') || isempty(filter_mode)
    is_bilinear = 1;
else
    is_bilinear = strcmp(filter_mode, 'bilinear');
end

for y_ul = 1:eh
    if mod(y_ul, 10) == 0
        fprintf('%d/%d\n', y_ul, eh);
    end
    y_c = y_ul + hw;
    for x_ul = 1:ew
        x_c = x_ul + hw;

        cosv =  cos_angle(y_c, x_c);
        sinv = -sin_angle(y_c, x_c);
        patch_x = x_c + du_grid * cosv + dv_grid *  sinv;
        patch_y = y_c + du_grid * sinv + dv_grid * -cosv;
        
        % 61, 74
        % 83, 92
        % 36, 98
        %if x_ul == 61 && y_ul == 74
        %if x_ul == 83 && y_ul == 92
        %    stophere;
        %end
        if is_bilinear
            patches(y_ul, x_ul, :) = lookup_bilinear(I, patch_x, patch_y);
        else
            patches(y_ul, x_ul, :) = lookup_nearest(I, patch_x, patch_y);
        end
    end
end
