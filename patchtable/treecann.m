function [nnf, T_total, runtime] = treecann(A, B, patch_w, A_grid, B_grid)

addpath 'TreeCann';
addpath 'TreeCANN/C_code';
addpath 'TreeCANN/matlab_tools/ann_wrapper/';

T0=clock; [nnf_dist, nnf_X, nnf_Y, runtime] = run_TreeCANN(A, B, patch_w, A_grid, B_grid); T_total = etime(clock, T0);

bew=size(B,2)-patch_w+1;
beh=size(B,1)-patch_w+1;
aew=size(A,2)-patch_w+1;
aeh=size(A,1)-patch_w+1;
nnf=zeros(aeh, aew, 3);
nnf(:,:,1)=min(max(nnf_X(1:aeh, 1:aew), 1), bew)-1;
nnf(:,:,2)=min(max(nnf_Y(1:aeh, 1:aew), 1), beh)-1;
nnf(:,:,3)=nnf_dist(1:aeh, 1:aew)/(255*255);
