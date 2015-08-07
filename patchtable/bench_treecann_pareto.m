maxNumCompThreads(1);
patch_w = 8;
anames = cell(0,0); bnames = cell(0, 0);
anames{length(anames)+1} = '../patchmatch/one_vs_many/a/001.png';
bnames{length(bnames)+1} = '../patchmatch/one_vs_many/b/001.png';
anames{length(anames)+1} = '../patchmatch/one_vs_many/a/002.png';
bnames{length(bnames)+1} = '../patchmatch/one_vs_many/b/002.png';
anames{length(anames)+1} = '../patchmatch/one_vs_many/a/003.png';
bnames{length(bnames)+1} = '../patchmatch/one_vs_many/b/003.png';
anames{length(anames)+1} = '../patchmatch/one_vs_many/a/004.png';
bnames{length(bnames)+1} = '../patchmatch/one_vs_many/b/004.png';
anames{length(anames)+1} = '../patchmatch/one_vs_many/a/005.png';
bnames{length(bnames)+1} = '../patchmatch/one_vs_many/b/005.png';
anames{length(anames)+1} = '../patchmatch/one_vs_many/a/006.png';
bnames{length(bnames)+1} = '../patchmatch/one_vs_many/b/006.png';
anames{length(anames)+1} = '../patchmatch/one_vs_many/a/007.png';
bnames{length(bnames)+1} = '../patchmatch/one_vs_many/b/007.png';
anames{length(anames)+1} = '../patchmatch/one_vs_many/a/008.png';
bnames{length(bnames)+1} = '../patchmatch/one_vs_many/b/008.png';
anames{length(anames)+1} = '../patchmatch/one_vs_many/a/009.png';
bnames{length(bnames)+1} = '../patchmatch/one_vs_many/b/009.png';
for i=1:length(anames)
    aimg{i} = imread(anames{i});
    bimg{i} = imread(bnames{i});
end
addpath 'TreeCANN';
addpath 'TreeCANN/C_code';
addpath 'TreeCANN/matlab_tools/ann_wrapper/';

itersL = [
1 2;
1 3;
1 4;
2 1;
2 3;
2 4;
2 5;
3 2;
2 7;
3 3;
3 4;
3 5;
3 6;
4 3;
3 8;
3 9;
4 4;
4 6;
4 7;
5 4;
4 10;
5 5;
5 6;
5 8;
6 6;
6 10;
8 7;
7 20;
9 15;
10 80];

%size(itersL)

for threads=[1]
    for idx=1:size(itersL,1)
        iters = itersL(idx,1);
        bgrid = itersL(idx,2);

        dL = []; tL = [];
        for i=1:length(anames)
            a = aimg{i};
            b = bimg{i};
            A=a; B=b;
            T0=clock; [nnf_dist, nnf_X , nnf_Y, runtime] =run_TreeCANN(uint8(A),uint8(B),patch_w,iters,bgrid, [], [], [], [], [], [], [], 1); T = etime(clock, T0);
        T = runtime(7);
            bew=size(B,2)-patch_w+1;
            beh=size(B,1)-patch_w+1;
            aew=size(A,2)-patch_w+1;
            aeh=size(A,1)-patch_w+1;
            nnf=zeros(aeh, aew, 3);
            nnf(:,:,1)=min(max(nnf_X(1:aeh, 1:aew), 1), bew)-1;
            nnf(:,:,2)=min(max(nnf_Y(1:aeh, 1:aew), 1), beh)-1;
            delete 'treecann_nnf.mat';
            save('treecann_nnf.mat', 'nnf');
            [status, out] = system(sprintf('/usr/bin/python check_dist.py %s %s treecann_nnf.mat', anames{i}, bnames{i}));
            davg = str2num(out);
            dL = [dL davg];
            tL = [tL T];
        end
        fprintf('%f %f "TREECANN(agrid=%d, bgrid=%d, threads=%d)"\n', mean(tL), mean(dL), iters, bgrid, threads);
    end
end
