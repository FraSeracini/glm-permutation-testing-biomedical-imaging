function t_values = Two_Sample_T_Stat(FA_group1, FA_group2, X, wm_mask, dim1, dim2, dim3)
%
% This function computes the voxel-wise two-sample t-statistic for comparing 
% FA (Fractional Anisotropy) values between two subject groups using the 
% General Linear Model (GLM).
%
% The GLM used is: Y = X1β1 + X2β2 + e
%
% INPUTS:
%   - FA_group1: 4D array (dim1 x dim2 x dim3 x n1) containing FA values for group 1
%   - FA_group2: 4D array (dim1 x dim2 x dim3 x n2) containing FA values for group 2
%   - X: Design matrix of size (n1 + n2) x 2, encoding group membership
%   - dim1, dim2, dim3: Dimensions of the 3D brain volume
%
% OUTPUT:
%   - t_values: 3D array (dim1 x dim2 x dim3) containing computed t-statistics for each voxel

% Initialize
t_values = zeros(dim1, dim2 , dim3);

% Compute the perpendicular projection operator
PX = X * pinv(X' * X) * X';

% Compute Rx
RX = eye(size(PX)) - PX;

% Compute the two-sample t-statistic 
for i = 1 : dim1
    for j = 1 : dim2
        for k = 1 : dim3
            if wm_mask(i, j, k) > 0

                % Total response
                Y = [squeeze(FA_group1(i, j, k, :)); squeeze(FA_group2(i, j, k, :))];

                % Determine e_hat
                e_hat = RX * Y;

                % Determine the estimate to the model parameters of the GLM
                beta_hat = pinv(X' * X) * X' * Y;
                
                % Estimate the variance of the stochastic component e_hat
                sigma2 = sum(e_hat.^2) / (size(Y,1) - rank(X));

                % Estimate the covariance matrix of beta_hat
                cov_beta = sigma2 * pinv(X' * X);

                % Contrast vector
                lambda = [1; -1];

                % Determine the t-statistic
                t_values(i, j, k) = (lambda' * beta_hat) ./ sqrt(lambda' * cov_beta * lambda);
            end
        end
    end
end
end