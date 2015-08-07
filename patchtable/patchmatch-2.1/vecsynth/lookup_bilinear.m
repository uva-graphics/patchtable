function out = lookup_bilinear(I, x, y)

I_size = size(I);
h = size(I, 1);
w = size(I, 2);

eps = 0.00001;
x = min(max(x, 1), w-eps);
y = min(max(y, 1), h-eps);

xi = floor(x);
yi = floor(y);
xf = x-xi;
yf = y-yi;

c_ul = I(sub2ind(I_size, yi,   xi  ));
c_ur = I(sub2ind(I_size, yi,   xi+1));
c_ll = I(sub2ind(I_size, yi+1, xi  ));
c_lr = I(sub2ind(I_size, yi+1, xi+1));

c_u = c_ul + (c_ur - c_ul) .* xf;
c_l = c_ll + (c_lr - c_ll) .* xf;

out = c_u + (c_l - c_u) .* yf;
