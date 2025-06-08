function [labels, mu, sigma2, mixtureProb, membershipProbMat] = fitGMM_1d_template(x, K)
% fitGMM_1d implements 1-dim Gaussian Mixture Modeling.

% Iteration params
maxIter = 500;
logL = nan(1, maxIter);
emtol = 1e-9; 

% Define objects
x = x(:);
N = length(x);
labels = nan(N, 1);
mu = nan(1, K);
sigma2 = nan(1, K);
mixtureProb = nan(1, K); 
l = nan(N, 1);
membershipProbMat = zeros(N, K);

% Set random initial values
rng(42);
mu = x(randi(N, 1, K))';
sigma2 = var(x) * ones(1, K); 
mixtureProb = ones(1, K) / K;


% EM iteration
for iter = 1:maxIter
    % E-step
    for k = 1:K
        membershipProbMat(:, k) = mixtureProb(k) * normpdf(x, mu(k), sqrt(sigma2(k)));
    end
    
    l = sum(membershipProbMat, 2);
    membershipProbMat = membershipProbMat ./ l; % normalize each sample
    
 
   

    % M-step
    Nrowsum = sum(membershipProbMat, 1);  % (1,3)

    mu = (membershipProbMat' * x) ./ Nrowsum'; %  (3,1500) * (1500, 1) / (3,1)
    % disp(size(mu)) % (3,1)
    sigma2 = (membershipProbMat' * (x.^2) - 2 .* mu' .* (membershipProbMat' * x) + mu'.^2 .* Nrowsum') ./ Nrowsum';
    mixtureProb = Nrowsum / N;
    
  

    % Terminate if converged
    if iter > 1
        deltaLogL = abs((logL(iter) - logL(iter-1)) /  logL(iter-1));
        if deltaLogL < emtol; disp(iter); disp(logL(iter)); break; end
    end
end

% Determine memberships
[~, labels] = max(membershipProbMat, [], 2);

end
