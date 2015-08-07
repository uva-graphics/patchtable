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
addpath 'CSH_code_v2';
cd CSH_code_v2; AddPaths; cd ..;
for threads=[1]
    for iters=30:-1:1
        dL = []; tL = [];
        for i=1:length(anames)
            a = aimg{i};
            b = bimg{i};
            T0=clock; nnf=CSH_nn(a,b,8,iters); T = etime(clock, T0);
            nnf = nnf(1:end-7,1:end-7,:);
            nnf1 = nnf-1;
            delete 'csh_nnf.mat';
            save('csh_nnf.mat', 'nnf1');
            [status, out] = system(sprintf('/usr/bin/python check_dist.py %s %s csh_nnf.mat', anames{i}, bnames{i}));
            davg = str2num(out);
            dL = [dL davg];
            tL = [tL T];
        end
        fprintf('%f %f "CSH(iters=%d, threads=%d)"\n', mean(tL), mean(dL), iters, threads);
    end
end
