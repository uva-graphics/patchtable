addpath ../patchmatch;

patch_w = 8;
vis = 0;
scale = 255;    % Use 255 for int mode, 1 for float/double mode

filename0 = 'lighthouse_small.png';
%filename0 = 'marriage_small.png';
I=im2double(imread(filename0));

close all;

for i=1:patch_w^2
    filename = 'out.pfm';
    if i > 1
        filename = sprintf('out%d.pfm', i);
    end
    a=readpfm(filename)/scale;

    
    K=load(sprintf('gck_%d_%03d.txt', patch_w, i));
    Ip=imfilter(I,fliplr(flipud(K)),'full',0);
    c = patch_w;
    if i == 1
        c = patch_w^2;
    end
    if vis
        figure;imshow(Ip/c)
        figure;imshow(a/c)
    end
    d=Ip-a;
    d=sqrt(sum(d.^2,3));
    if vis
        figure;imagesc(d);
    end
    fprintf('Kernel %d difference: %f\n', i, max(d(:)));
end
