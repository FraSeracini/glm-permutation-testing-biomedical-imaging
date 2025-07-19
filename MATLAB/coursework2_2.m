%% Part 2
%% 2
%% a


% Load ROI mask
fid = fopen('wm_mask.img', 'r', 'l');
wm_mask = fread(fid, 'float');
fclose(fid);
wm_mask = reshape(wm_mask, [40, 40, 40]);

% Number of subjects
n1 = 8;  
n2 = 8;  
n = n1 + n2;

% Load FA images
FA_group1 = zeros(40,40,40,n1);
FA_group2 = zeros(40,40,40,n2);

for i = 1:n1
    filename = sprintf('CPA%d_diffeo_fa.img', i); 
    fid = fopen(filename, 'r', 'l');
    data = fread(fid, 'float');
    fclose(fid);
    FA_group1(:,:,:,i) = reshape(data, [40, 40, 40]);
end

for i = 1:n2
    filename = sprintf('PPA%d_diffeo_fa.img', i);
    fid = fopen(filename, 'r', 'l');
    data = fread(fid, 'float');
    fclose(fid);
    FA_group2(:,:,:,i) = reshape(data, [40, 40, 40]);
end

t_values = zeros(40, 40 ,40);

X = [ones(n1,1), zeros(n1,1); zeros(n2,1), ones(n2,1)];
% Compute projection matrix PX
PX = X * pinv(X' * X) * X';
RX = eye(size(PX)) - PX;

for i = 1 : 40
    for j = 1 : 40
        for k = 1 : 40
            if wm_mask(i, j, k) > 0
                Y = [squeeze(FA_group1(i, j, k, :)); squeeze(FA_group2(i, j, k, :))];
                e_hat = RX * Y;
                % Compute least squares estimate of beta
                beta_hat = pinv(X' * X) * X' * Y;
                % Estimate MSE
                sigma_sq = sum(e_hat.^2) / (size(Y,1) - rank(X));

                % Compute covariance matrix of beta_hat
                S_beta_hat = sigma_sq * pinv(X' * X);

                % Compute contrast vector
                lambda = [1; -1];
                t_values(i, j, k) = (lambda' * beta_hat) ./ sqrt(lambda' * S_beta_hat * lambda);
            end
        end
    end
end

% Compute original max t-statistic
[max_val, linear_idx] = max(t_values(:));
[i_max, j_max, k_max] = ind2sub(size(t_values), linear_idx);
%% b

perms = nchoosek(1:n, n1);
num_permutations = size(perms, 1);
max_t_val = zeros(num_permutations, 1);
D = cat(4, FA_group1, FA_group2);


for p = 1:num_permutations
    perm_idx = perms(p, :);
    FA_group1 = D(:, :, :, perm_idx);
    FA_group2 = D(:, :, :, setdiff(1:n, perm_idx));
    for i = 1 : 40
        for j = 1 : 40
            for k = 1 : 40
                if wm_mask(i, j, k) > 0
                    Y = [squeeze(FA_group1(i, j, k, :)); squeeze(FA_group2(i, j, k, :))];
                    e_hat = RX * Y;
                    % Compute least squares estimate of beta
                    beta_hat = pinv(X' * X) * X' * Y;
                    % Estimate MSE
                    sigma_sq = sum(e_hat.^2) / (size(Y,1) - rank(X));

                    % Compute covariance matrix of beta_hat
                    S_beta_hat = sigma_sq * pinv(X' * X);

                    % Compute contrast vector
                    lambda = [1; -1];
                    t_values(i, j, k) = (lambda' * beta_hat) ./ sqrt(lambda' * S_beta_hat * lambda);
                end
            end
        end
    end
    max_t_val(p) = max(t_values(:));
end
%% 

histogram(max_t_val, 50); % 50 bin per visualizzare la distribuzione
xlabel('Maximum t-statistic');
ylabel('Frequency');
title('Empirical Distribution of Maximum t-statistic');

p_value_corrected = sum(max_t_val >= max_val) / length(max_t_val)

threshold_95 = prctile(max_t_val, 95)



