%% Econometric Methods I Problem Set 4
% Replicating Alsan: 2. Linear regression with controls 
% Barcelona School of Economics, 2025-2026
% Author: Lea Röller

%% housekeeping
clear all; close all; clc; 

cd() %specify folder path if necessary

%% load data set
load(fullfile('..','00 Clean data','precolonial.mat'));

clearvars -except animals intensive plow female_ag ln_popd_murdock slavery_indigenous central tse ...
                  meanrh meantemp itx prop_tropics malaria_index river coast lon abslat meanalt SI  %keep variables that we want and drop the rest for now


%% Question one: linear regression

% List of dependent variables
depVars = {'animals', 'intensive', 'plow', 'female_ag', ...
           'ln_popd_murdock', 'slavery_indigenous', 'central'};

% adding the controls!
controls = {'meanrh', 'meantemp', 'itx', 'prop_tropics', ...
            'malaria_index', 'river', 'coast', 'lon', 'abslat', 'meanalt', 'SI'};

% Initialize result table
results = table('Size', [numel(depVars), 8], ...
    'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'string', 'double'}, ...
    'VariableNames', {'Variable', 'b2', 's2', 'SE_b2', 't_stat', 'p_value', 'Decision', 'n_obs'});

alfa = 0.05;


for i = 1:numel(depVars)
    y = eval(depVars{i});
    
    % Extract tse and all controls
    x = tse;  % main independent variable
    for c = 1:numel(controls)
        x = [x, eval(controls{c})];  % append control variable as additional column
    end

    % Remove rows with missing y or any missing controls
    valid = ~isnan(y) & all(~isnan(x), 2);
    y = y(valid);
    x = x(valid, :);

    X = [ones(length(y), 1), x]; % Construct regression matrix with intercept
    beta_hat = (X' * X) \ (X' * y);  % OLS estimation

    % Residuals and variance
    n = length(y); % we need to take the actual length of y that we have for the SEs!
    k = size(X, 2);
    residuals = y - X * beta_hat;
    sigma2_hat = (residuals' * residuals) / (n - k);
    cov_beta = sigma2_hat * inv(X' * X);
    se_beta = sqrt(diag(cov_beta));

    % t-test for tse coefficient
    t_beta2 = beta_hat(2) / se_beta(2);
    p_beta2 = 2 * (1 - tcdf(abs(t_beta2), n - k)); % two-sided

    % Decision
    decision2 = "Fail to reject $H_0$";
    if p_beta2 < alfa
        decision2 = "Reject $H_0$";
    end

    % Store results
    results.Variable(i) = depVars{i};
    results.n_obs(i) = n;
    results.b2(i) = beta_hat(2);
    results.s2(i) = sigma2_hat;
    results.SE_b2(i) = se_beta(2);
    results.t_stat(i) = t_beta2;
    results.p_value(i) = p_beta2;
    results.Decision(i) = decision2;
end

% Display table
disp(results);

%% Export the results
% File name to save
% File name to save
outputFolder = fullfile('..','02 Outputs');
filename = fullfile(outputFolder, 'Exercise_1.2.tex');  % full path
fid = fopen(filename, 'w'); % Open file for writing

% Write LaTeX table header
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{lrrrrrrrl}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'Variable & $b_2$ & $SE(b_2)$ & t-stat & p-value & Decision & $s^2$ & Observations \\\\\n');
fprintf(fid, '\\hline\n');

% Loop through results and write each row
% Loop through results and write each row with variable names in $...$
for i = 1:height(results)
    fprintf(fid, '$%s$ & %.3f & %.3f & %.3f & %.3f & %s & %.3f & %d\\\\\n', ...
        results.Variable(i), ...
        results.b2(i), ...
        results.SE_b2(i), ...
        results.t_stat(i), ...
        results.p_value(i), ...
        results.Decision(i), ...
        results.s2(i), ...
        results.n_obs(i));  % only keep this one for Observations
end

% Write LaTeX table footer
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{OLS regression results for all outcome variables}\n');
fprintf(fid, '\\label{tab:regression_results}\n');

% Close the file
fclose(fid);


%% Cross-check

% Example for one dependent variable
y = central;   % choose which variable to check

% List of control variables
controls = {'meanrh', 'meantemp', 'itx', 'prop_tropics', ...
            'malaria_index', 'river', 'coast', 'lon', 'abslat', 'meanalt', 'SI'};

% Extract tse and controls
x = tse;  % main independent variable
for c = 1:numel(controls)
    x = [x, eval(controls{c})];  % append control variables
end

% Remove rows with missing y or any control
valid = ~isnan(y) & all(~isnan(x), 2);
y = y(valid);
x = x(valid, :);

% Prepare table for fitlm
varNames = ['Y', 'TSE', controls];   % variable names for table
tbl = array2table([y, x], 'VariableNames', varNames);

% Construct formula: Y ~ TSE + control1 + control2 + ...
formula = ['Y ~ TSE'];
for i = 1:numel(controls)
    formula = [formula, ' + ', controls{i}];
end

% Run built-in linear regression
mdl = fitlm(tbl, formula);

% Display regression table
disp(mdl)

%%