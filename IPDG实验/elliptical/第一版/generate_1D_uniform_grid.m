function  [h,mid_points]=generate_1D_uniform_grid(left,right,ng)

h = (right-left)/ng;
mid_points = linspace(left+h/2,right-h/2,ng);
h = ones(1,ng)*h;
