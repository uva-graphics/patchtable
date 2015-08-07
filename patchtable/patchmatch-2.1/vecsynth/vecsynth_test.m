
% Exemplar images
itest = 1:1;
target_size = 128; %32; %128; %128; %32; %128;
addpath ..;

img = cell(9, 1);
for i=1:length(img)
    img{i} = im2double(imread(sprintf('images/sdf_%d.png', i-1)));
    img{i} = rgb2gray(imresize(img{i}, [target_size target_size], 'bilinear'));
end

% Query img{i} versus db{i}
db = cell(length(img), 1);
for i=itest
    for j=1:length(img)
        if j==i
            continue
        end
        if isempty(db{i})
            db{i} = img{j};
        else
            db{i} = [db{i} img{j}];
        end
    end
end

% --------------------------------------------------------------------------------
% Rotations search with generalized PatchMatch
% --------------------------------------------------------------------------------

% Uses nearest neighbor filtering
% Exemplar patches are rotated; output (synthesized) patch is always upright

patch_w = 5;
cores = 8;
%nn_iters = 175;    % Nearest neighbor filtering
nn_iters = 28;
rs_max = [];
rs_min = [];
rs_ratio = [];
rs_iters = [];
bmask = [];
win_size = [];
nnfield_prev = [];
nnfield_prior = [];
prior_winsize = [];
knn = [];
scalerange = 1;

tic;

avg_err_gpm = [];
avg_time_gpm = [];
for i=itest
    T0 = toc;
    ann = nnmex(cat(3, img{i}, img{i}, img{i}), cat(3, db{i}, db{i}, db{i}), 'rotscale', patch_w, ...
                nn_iters, rs_max, rs_min, rs_ratio, rs_iters, cores, bmask, ...
                win_size, nnfield_prev, nnfield_prior, prior_winsize, knn, scalerange);
    T1 = toc;
    err_channel = sqrt(ann(1:end-patch_w+1,1:end-patch_w+1,3));
    avg_time_gpm = [avg_time_gpm (T1-T0)];
    avg_err_gpm = [avg_err_gpm mean(err_channel(:))];
end

avg_err_gpm
avg_time_gpm

% --------------------------------------------------------------------------------
% Rotations search with principal orientation + PCA
% --------------------------------------------------------------------------------

filter_mode = 'bilinear';
num_pca = 8;
avg_time_princomp = [];
avg_err_princomp = [];
nn_iters = 2;

view_pca = 1;

for i=itest
    [~,img_orient{i}] = imgradient(img{i});
    [~,db_orient{i}] = imgradient(db{i});
   
    img_patches{i} = oriented_patches(img{i}, img_orient{i}, patch_w, filter_mode);
    db_patches{i} = oriented_patches(db{i}, db_orient{i}, patch_w, filter_mode);
   
    img_patches_flat = reshape(img_patches{i}, [size(img_patches{i}, 1)*size(img_patches{i},2), size(img_patches{i}, 3)]);
    db_patches_flat  = reshape(db_patches{i},  [size(db_patches{i}, 1) *size(db_patches{i},2),  size(db_patches{i}, 3)]);

    all_patches_flat = [img_patches_flat; db_patches_flat];
    assert (size(all_patches_flat, 2) == patch_w^2);
    all_patches_flat = all_patches_flat - repmat(mean(all_patches_flat), [size(all_patches_flat, 1) 1]);   % Center
    [coeff, score, ~, ~, explained] = pca(all_patches_flat);
    fprintf('PCA variance explained (full list):\n');
    explained
    fprintf('Variance explained by first %d: %f\n', num_pca, sum(explained(1:num_pca)));
    %'stop 1';
    coeff = coeff(:,1:num_pca);    % score * coeff' reconstructs patches.    coeff' size: num_pca x patch_w^2
    score = score(:,1:num_pca);
    pca_reduced = score * coeff';

    img_patches_pca{i} = pca_reduced(1:size(img_patches_flat,1), 1:patch_w^2);
    db_patches_pca{i} = pca_reduced(size(img_patches_flat,1)+1:end, 1:patch_w^2);
    img_patches_pca{i} = reshape(img_patches_pca{i}, size(img_patches{i}));
    db_patches_pca{i} = reshape(db_patches_pca{i}, size(db_patches{i}));

    ew = size(img{i}, 2)-patch_w+1;
    eh = size(img{i}, 1)-patch_w+1;
    if view_pca
        figure;
        title('PCA bases');
        for j=1:num_pca
            subplot(1, num_pca, j);
            imagesc(reshape(coeff(:, j), [patch_w patch_w])');
            title(sprintf('Dim %d', j));
        end
    end
    
    T0 = toc;
    ann = nnmex(img_patches_pca{i}, db_patches_pca{i}, 'cputiled', 1, ...
                nn_iters, rs_max, rs_min, rs_ratio, rs_iters, cores, bmask, ...
                win_size, nnfield_prev, nnfield_prior, prior_winsize, knn, scalerange);
    T1 = toc;
    
    err_channel = zeros(ew, eh);
    for y=1:eh
        for x=1:ew
            xp = ann(y,x,1)+1;
            yp = ann(y,x,2)+1;
            delta = img_patches{i}(y,x,:) - db_patches{i}(yp,xp,:);
            err_channel(y,x) = sqrt(dot(delta, delta));
        end
    end
    %err_channel = sqrt(ann(:,:,3));
    avg_time_princomp = [avg_time_princomp T1-T0];
    avg_err_princomp = [avg_err_princomp mean(err_channel(:))*256];
end

avg_time_gpm
avg_err_gpm

avg_time_princomp
avg_err_princomp

'stop';