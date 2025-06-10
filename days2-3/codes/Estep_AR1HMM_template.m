function [c, gammaMat, xiArr] = ...
            Estep_AR1HMM_template(x, x1, initPr, tranPr, phi0, phi1, sigmasq, T, M)
% Estep_AR1HMM
%
% J Noh, 2025/02

%% Define objects
alphaMat = zeros(T, M);
betaMat = zeros(T, M);
gammaMat = zeros(T, M);
xiArr = zeros(T-1, M, M);
bMat = zeros(T, M); 

c = zeros(T, 1);              % normalizing scale factor for alpha_t(i)
d = zeros(T, 1);              % normalizing scale factor for beta_t(i)

%% pdf function values, b_i(x_t | x_(t-1))
for i = 1:M
    mu = phi0(i) + phi1(i) * x1(:);
    sigma = ( sigmasq(i) )^0.5;
    bMat(:, i) = pdf('Normal', x(:), mu(:), sigma);
end

%% alpha_t(i) forward equation
alphaMat(1, :) = initPr(:)' .* bMat(1, :);
c(1) = sum(alphaMat(1, :));
alphaMat(1, :) = alphaMat(1, :) / c(1); 

for t = 2:T
    for j = 1:M
        alphaMat(t, j) = sum(alphaMat(t-1, :) .* tranPr(:, j)') * bMat(t, j);
    end
    c(t) = sum(alphaMat(t, :));
    alphaMat(t, :) = alphaMat(t, :) / c(t);  % scaling
end


%% beta_t(j) backward equation
betaMat(T, :) = 1;          % unscaled initialization
d(T) = sum(betaMat(T, :));  % should be 1, but included for completeness
betaMat(T, :) = betaMat(T, :) / d(T);

for t = T-1:-1:1
    for i = 1:M
        betaMat(t, i) = sum(tranPr(i, :) .* bMat(t+1, :) .* betaMat(t+1, :));
    end
    d(t) = sum(betaMat(t, :));            % separate normalization factor
    betaMat(t, :) = betaMat(t, :) / d(t); % normalize
end


%% define gamma_t(i) 

for t = 1:T
    denom = sum(alphaMat(t, :) .* betaMat(t, :));
    for i = 1:M
        gammaMat(t, i) = (alphaMat(t, i) * betaMat(t, i)) / denom;
    end
end


%% define xi_t(i,j)  
for t = 1:T-1
    x_ = 0;
    for i = 1:M
        for j = 1:M
            x_ = x_ + alphaMat(t, i) * tranPr(i, j) * bMat(t+1, j) * betaMat(t+1, j);
        end
    end
    for i = 1:M
        for j = 1:M
            xiArr(t, i, j) = (alphaMat(t, i) * tranPr(i, j) * bMat(t+1, j) * betaMat(t+1, j)) / x_;
        end
    end



end
