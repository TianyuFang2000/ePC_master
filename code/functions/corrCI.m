function [r_lower, r_upper] = corrCI(r, n, alpha)
    % corrCI computes the confidence interval for a Pearson correlation coefficient
    % Inputs:
    %   r     - observed Pearson correlation coefficient
    %   n     - number of paired observations
    %   alpha - significance level (e.g., 0.05 for 95% CI)
    % Outputs:
    %   r_lower, r_upper - lower and upper bounds of the confidence interval

    if nargin < 3
        alpha = 0.05;  % default to 95% CI
    end

    % Fisher r-to-z transformation
    z = 0.5 * log((1 + r) / (1 - r));
    
    % Standard error of z
    SE = 1 / sqrt(n - 3);
    
    % Z critical value (two-tailed)
    z_crit = norminv(1 - alpha / 2);
    
    % Confidence interval in z-space
    z_lower = z - z_crit * SE;
    z_upper = z + z_crit * SE;
    
    % Convert back to r-space
    r_lower = (exp(2*z_lower) - 1) / (exp(2*z_lower) + 1);
    r_upper = (exp(2*z_upper) - 1) / (exp(2*z_upper) + 1);
end
