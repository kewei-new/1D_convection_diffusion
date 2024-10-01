function xnd = nd(x)

% Defining constants
ya = 5e5;
yb = 2e3;

xjul = 0.1;
xjur = 0.5;
xjwd = 0.06;
xjwh = xjwd / 2;
xjull = xjul - xjwh;
xjulr = xjul + xjwh;
xjurl = xjur - xjwh;
xjurr = xjur + xjwh;

% Calculating xnd based on x00
if x < xjull
    xnd = ya;
elseif x < xjulr
    yr = (x - xjull) / (xjwd + 1e-20);
    xnd = (ya - yb) * (1 - yr^3)^3 + yb;
elseif x < xjurl
    xnd = yb;
elseif x < xjurr
    yr = (x - xjurl) / (xjwd + 1e-20);
    xnd = (yb - ya) * (1 - yr^3)^3 + ya;
else
    xnd = ya;
end