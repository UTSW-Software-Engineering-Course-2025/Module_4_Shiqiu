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

% Set random initial values




% EM iteration
for iter = 1:maxIter
    % M-step
    
    
    
    

    % E-step
    
    
    
    
    
    % Terminate if converged
    if iter > 1
        deltaLogL = abs((logL(iter) - logL(iter-1)) /  logL(iter-1));
        if deltaLogL < emtol; disp(iter); disp(logL(iter)); break; end
    end
end

% Determine memberships








end
