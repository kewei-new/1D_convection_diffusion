function doping = nd_profile(x, params)
%ND_PROFILE Smooth p-n doping profile used by the device mode.

    dop = params.device.doping;

    ya = dop.ya;
    yb = dop.yb;
    x_l = dop.x_l;
    x_r = dop.x_r;
    xjwd = dop.transition_width;
    xjwh = xjwd / 2;

    x_ll = x_l - xjwh;
    x_lr = x_l + xjwh;
    x_rl = x_r - xjwh;
    x_rr = x_r + xjwh;

    doping = zeros(size(x));

    mask1 = (x < x_ll);
    mask2 = (x >= x_ll) & (x < x_lr);
    mask3 = (x >= x_lr) & (x < x_rl);
    mask4 = (x >= x_rl) & (x < x_rr);
    mask5 = (x >= x_rr);

    doping(mask1) = ya;

    yr = (x(mask2) - x_ll) ./ (xjwd + 1e-20);
    doping(mask2) = (ya - yb) .* (1 - yr.^3).^3 + yb;

    doping(mask3) = yb;

    yr = (x(mask4) - x_rl) ./ (xjwd + 1e-20);
    doping(mask4) = (yb - ya) .* (1 - yr.^3).^3 + ya;

    doping(mask5) = ya;
end
