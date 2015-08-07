function out = lookup_nearest(I, x, y)

I_size = size(I);
h = size(I, 1);
w = size(I, 2);

xi = floor(x+0.5);
yi = floor(y+0.5);

xi = min(max(xi, 1), w);
yi = min(max(yi, 1), h);

out = I(sub2ind(I_size, yi, xi));
