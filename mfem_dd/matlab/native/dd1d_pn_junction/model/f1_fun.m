function result = f1_fun(X,t)

result = -sin(t).*sin(X)+cos(t).*sin(X)-2*sin(t).*cos(t).*sin(X).*cos(X);

% result = -sin(t).*sin(X)+cos(t).*sin(X);

% result = 0*X + 0*t;