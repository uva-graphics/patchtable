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
    aimg{i} = im2double(imread(anames{i}));
    bimg{i} = im2double(imread(bnames{i}));
end
addpath 'patchmatch-2.1';
for threads=[1]
    for iters=[[1 2 3 4 5 6 7 8 9 10] 15:5:100]
        dL = []; tL = [];
        for i=1:length(anames)
            a = aimg{i};
            b = bimg{i};
            tic; nnf=nnmex(a,b,'cputiled',8,iters,[],[],[],[],threads); T = toc;
            %d=nnf(:,:,3); dvalid=d(1:end-8,1:end-8,:); davg = mean(sqrt(double(dvalid(:))/255^2));
            nnf_d = double(nnf);
            nnf_d = nnf_d(1:end-7,1:end-7,:);
            save('pm_nnf.mat', 'nnf_d');
            [status, out] = system(sprintf('/usr/bin/python check_dist.py %s %s pm_nnf.mat', anames{i}, bnames{i}));
            davg = str2num(out);
            dL = [dL davg];
            tL = [tL T];
        end
        fprintf('%f %f "PM(iters=%d, threads=%d)"\n', mean(tL), mean(dL), iters, threads);
    end
end
