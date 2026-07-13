function Gauss_coefficient = generate_1D_gaussian_quadrature_r(n)
%GENERATE_1D_GAUSSIAN_QUADRATURE_R
% Gauss-Legendre quadrature on reference interval [-1/2, 1/2].
%
% Output:
%   Gauss_coefficient(:,1) : quadrature points
%   Gauss_coefficient(:,2) : quadrature weights
%
% This function is added as a compatibility wrapper for old code.

    switch n
        case 1
            x = 0;
            w = 1;

        case 2
            x_std = [-1/sqrt(3); 1/sqrt(3)];
            w_std = [1; 1];

            x = x_std / 2;
            w = w_std / 2;

        case 3
            x_std = [-sqrt(3/5); 0; sqrt(3/5)];
            w_std = [5/9; 8/9; 5/9];

            x = x_std / 2;
            w = w_std / 2;

        case 4
            x_std = [-0.8611363115940526;
                     -0.3399810435848563;
                      0.3399810435848563;
                      0.8611363115940526];
            w_std = [0.3478548451374538;
                     0.6521451548625461;
                     0.6521451548625461;
                     0.3478548451374538];

            x = x_std / 2;
            w = w_std / 2;

        case 5
            x_std = [-0.9061798459386640;
                     -0.5384693101056831;
                      0;
                      0.5384693101056831;
                      0.9061798459386640];
            w_std = [0.2369268850561891;
                     0.4786286704993665;
                     0.5688888888888889;
                     0.4786286704993665;
                     0.2369268850561891];

            x = x_std / 2;
            w = w_std / 2;

        otherwise
            error('generate_1D_gaussian_quadrature_r: unsupported n = %d', n);
    end

    Gauss_coefficient = [x, w];
end