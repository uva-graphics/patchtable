addpath /Users/connelly/filter-approx/proj/patchmatch

ann=readpfm('/Users/connelly/filter-approx/proj/patchdb/out.pfm');
a=im2double(imread('a.png'));
b=im2double(imread('b.png'));
ann_pad=int32(zeros(size(a,1),size(a,2),3));
ann_pad(1:size(ann,1),1:size(ann,2),:) = ann;

cores = 1;

tic;nnf=nnmex(a,b,'cputiled',8,1,0,[],[],[],cores,[],[],int32(ann_pad));toc
d=nnf(:,:,3); dvalid=d(1:end-8,1:end-8,:); mean(sqrt(double(dvalid(:))/255^2))
