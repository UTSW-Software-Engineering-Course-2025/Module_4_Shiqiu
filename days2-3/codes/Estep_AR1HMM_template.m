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



%% beta_t(j) backward equation



%% define gamma_t(i) 



%% define xi_t(i,j)  



end
