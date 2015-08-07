
% Simple visualization script that runs PatchTable and compares also with PatchMatch

addpath ../patchmatch
fprintf('PatchTable query time is reported below as "1-NN overall lookup time."\n\n');
if ispc
    patchtable = 'patchtable';
else
    patchtable = './patchtable';
end
system(sprintf('%s match vidpair0/a.png vidpair0/b.png out.pfm -ndims 6 -limit 340 -do_rs 1 -parallel_dt 1 -randomize_dt 1 -dt_threads 8 -threads 1', patchtable));
nnf=readpfm('out.pfm');
close all;
figure('name', 'PatchTable NNF x coord'); imagesc(nnf(:,:,1))
figure('name', 'PatchTable NNF y coord'); imagesc(nnf(:,:,2))
figure('name', 'PatchTable NNF dist'); imagesc(sqrt(nnf(:,:,3)))

addpath patchmatch-2.1
a=im2double(imread('vidpair0/a.png'));
b=im2double(imread('vidpair0/b.png'));
fprintf('PatchMatch time:\n');
tic; nnf=nnmex(a,b,'cputiled',8,5,[],[],[],[],1);toc;
figure('name', 'PatchMatch NNF x coord');
imagesc(nnf(:,:,1))
fprintf('PatchMatch mean dist:\n');
d=nnf(:,:,3); dvalid=d(1:end-8,1:end-8,:); mean(sqrt(double(dvalid(:))/255^2))
figure('name', 'PatchMatch NNF y coord');
imagesc(nnf(:,:,2))
figure('name', 'PatchMatch NNF dist');
imagesc(sqrt(double(dvalid)))

