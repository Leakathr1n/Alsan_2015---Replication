%% Econometric Methods I Problem Set 4
% Replicating Alsan: 3. The differential effect of the TseTse Suitability Index in tropical Africa  
% Barcelona School of Economics, 2025-2026
% Author: Lea Röller

%% housekeeping
clear all; close all; clc; 

cd() %specify folder path if necessary

%% load the data
load(fullfile('..','00 Clean data','placebo.mat'));

% List of variables to always keep
keepVars = {'animals', 'intensive', 'plow', 'female_ag', ...
            'slavery_indigenous', 'central', 'tse', ...
            'meanrh', 'meantemp', 'itx', 'malaria_index', ...
            'river', 'coast', 'lon', 'abslat', 'meanalt', 'SI'};

% Find all variables starting with 'africa_'
allVars = who; % list all variables in workspace
africaVars = allVars(startsWith(allVars, 'africa'));

keepVars = [keepVars, africaVars']; % Combine the lists
clearvars('-except', keepVars{:}); % Clear all other variables

%% Question one: regression

% List of dependent variables
depVars = {'animals', 'intensive', 'plow', 'female_ag', ...
             'slavery_indigenous', 'central'};

% Controls
controls = {'meanrh', 'meantemp', 'itx', ...
            'malaria_index', 'river', 'coast', 'lon', 'abslat', 'meanalt', 'SI'};

% Initialize result table
results = table('Size', [numel(depVars), 12], ...
    'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'string', 'double'}, ...
    'VariableNames', {'Variable', 'b2_tse', 'SE_tse', 't_tse', 'p_tse', ...
                      'b2_africaTSE', 'SE_africaTSE', 't_africaTSE', 'p_africaTSE', ...
                      's2', 'Decision', 'n_obs'});

alfa = 0.05;

for i = 1:numel(depVars)
    y = eval(depVars{i});
    
    % Start with main variables
    x = [tse, africa];  % main variables

    % Include pre-computed interaction terms (africa_x)
    for c = 1:numel(controls)
        ctrl = eval(controls{c});           % control
        x = [x, ctrl];                       % control
        x = [x, eval(['africa_', controls{c}])]; % precomputed interaction control*africa
    end

    % Include africa*tse if precomputed
    x = [x, africa_tse]; 

    % Remove rows with missing y or any missing x
    valid = ~isnan(y) & all(~isnan(x), 2);
    y = y(valid);
    x = x(valid, :);

    % Regression matrix with intercept
    X = [ones(length(y), 1), x];
    beta_hat = (X' * X) \ (X' * y);

    % Residuals and variance
    n = length(y);
    k = size(X, 2);
    residuals = y - X * beta_hat;
    sigma2_hat = (residuals' * residuals) / (n - k);
    cov_beta = sigma2_hat * inv(X' * X);
    se_beta = sqrt(diag(cov_beta));

    % t-test for tse coefficient
    t_tse = beta_hat(2) / se_beta(2);
    p_tse = 2 * (1 - tcdf(abs(t_tse), n - k));

    % t-test for africa*tse coefficient (assuming column position is last)
    t_africaTSE = beta_hat(end) / se_beta(end);
    p_africaTSE = 2 * (1 - tcdf(abs(t_africaTSE), n - k));

    % Decision based on tse coefficient
    decision = "Fail to reject $H_0$";
    if p_tse < alfa
        decision = "Reject $H_0$";
    end

    % Decision for africa*tse coefficient
    decision_africaTSE = "Fail to reject $H_0$";
    if p_africaTSE < alfa
        decision_africaTSE = "Reject $H_0$";
    end

    % Store results
    results.Variable(i) = depVars{i};
    results.n_obs(i) = n;
    results.s2(i) = sigma2_hat;

    results.b2_tse(i) = beta_hat(2);
    results.SE_tse(i) = se_beta(2);
    results.t_tse(i) = t_tse;
    results.p_tse(i) = p_tse;

    results.b2_africaTSE(i) = beta_hat(end);
    results.SE_africaTSE(i) = se_beta(end);
    results.t_africaTSE(i) = t_africaTSE;
    results.p_africaTSE(i) = p_africaTSE;

    results.Decision(i) = decision;
    results.Decision_africaTSE(i) = decision_africaTSE;
end

disp(results); % Display table


%% Export to LaTeX

% File name for TSE coefficients
filenameTSE = fullfile(outputFolder, 'Exercise_1.3_TSE.tex');
fid = fopen(filenameTSE, 'w');

fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{lrrrrrl}\n'); % Variable + b, SE, t, p, s^2, Decision = 6 columns
fprintf(fid, '\\hline\n');
fprintf(fid, 'Variable & $\\beta_2$ & $SE(\\beta_2)$ & t-stat & p-value & s^2 & Decision \\\\\n');
fprintf(fid, '\\hline\n');

for i = 1:height(results)
    fprintf(fid, '$%s$ & %.2f & %.2f & %.2f & %.2f & %.2f & %s \\\\\n', ...
        results.Variable(i), ...
        results.b2_tse(i), ...
        results.SE_tse(i), ...
        results.t_tse(i), ...
        results.p_tse(i), ...
        results.s2(i), ...
        results.Decision(i));
end

fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{OLS regression results: TSE coefficient only, including controls.}\n');
fprintf(fid, '\\label{tab:regression_results_TSE}\n');

fclose(fid);

% File name for interaction coefficients
filenameInt = fullfile(outputFolder, 'Exercise_1.3_TSE_Africa.tex');
fid = fopen(filenameInt, 'w');

fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{lrrrrrl}\n'); % Variable + b, SE, t, p, s^2, Decision = 6 columns
fprintf(fid, '\\hline\n');
fprintf(fid, 'Variable & $\\gamma_2$ & $SE(\\gamma_2)$ & t-stat & p-value & s^2 & Decision \\\\\n');
fprintf(fid, '\\hline\n');

for i = 1:height(results)
    fprintf(fid, '$%s$ & %.2f & %.2f & %.2f & %.2f & %.2f & %s \\\\\n', ...
        results.Variable(i), ...
        results.b2_africaTSE(i), ...
        results.SE_africaTSE(i), ...
        results.t_africaTSE(i), ...
        results.p_africaTSE(i), ...
        results.s2(i), ...
        results.Decision_africaTSE(i));
end

fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{OLS regression results: TSE*Africa interaction coefficient only, including controls.}\n');
fprintf(fid, '\\label{tab:regression_results_TSE_Africa}\n');

fclose(fid);

disp('LaTeX table exported!');

%% Question 1.2 - an F-test - DOUBLE - CHECK THIS AGAIN!!

% List of dependent variables
depVars = {'animals', 'intensive', 'plow', 'female_ag', ...
             'slavery_indigenous', 'central'};

% Controls
controls = {'meanrh', 'meantemp', 'itx', ...
            'malaria_index', 'river', 'coast', 'lon', 'abslat', 'meanalt', 'SI'};

% Initialize table for linear combination tests
FtestResults = table('Size', [numel(depVars), 4], ...
    'VariableTypes', {'string','double','double','string'}, ...
    'VariableNames', {'Variable','F_stat','p_value','Decision'});

alfa = 0.05;

for i = 1:numel(depVars)
    y = eval(depVars{i});
    
    % Start with main variables
    x = [tse, africa];                 % main variables
    x = [x, africa_tse];               % precomputed TSE*Africa

    % Add controls and interactions with Africa
    for c = 1:numel(controls)
        ctrl = eval(controls{c});
        x = [x, ctrl];                               % control
        x = [x, eval(['africa_', controls{c}])];    % precomputed interaction
    end

    % Remove rows with missing y or any missing x
    valid = ~isnan(y) & all(~isnan(x), 2);
    y = y(valid);
    x = x(valid, :);

    % Regression matrix with intercept
    X = [ones(length(y), 1), x];
    beta_hat = (X' * X) \ (X' * y);

    % Residuals and covariance
    n = length(y);
    k = size(X, 2);
    residuals = y - X * beta_hat;
    sigma2_hat = (residuals' * residuals) / (n - k);
    cov_beta = sigma2_hat * inv(X' * X);

    % Linear combination test: beta2 + interaction = 0
    R = zeros(1, k);
    R(2) = 1;  % beta2 (TSE)
    R(4) = 1;  % beta4 (TSE*Africa)
    r = 0;

    df_num = rank(R);        % numerator df
    df_den = n - k;          % denominator df

    F = (R*beta_hat - r)' * inv(R * sigma2_hat * inv(X'*X) * R') * (R*beta_hat - r);
    p_value = 1 - fcdf(F, df_num, df_den);

    % Decision
    decisionF = "Fail to reject $H_0$"; % write in latex format so that we can easily export
    if p_value < alfa
        decisionF= "Reject $H_0$";
    end

    % Store results
    FtestResults.Variable(i) = depVars{i};
    FtestResults.F_stat(i) = F;
    FtestResults.p_value(i) = p_value;
    FtestResults.Decision(i) = decisionF;
end

% Display table
disp(FtestResults);

%% Export to LaTeX

% File name to save
outputFolder = fullfile('..','02 Outputs');
filename = fullfile(outputFolder, 'Exercise_1.3_Ftest.tex');  % full path
fid = fopen(filename, 'w'); % Open file for writing

% LaTeX table header
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{lrrl}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'Variable & F-stat & p-value & Decision \\\\\n');
fprintf(fid, '\\hline\n');

% Loop through results and write each row (2 decimals for numeric values)
for i = 1:height(FtestResults)
    fprintf(fid, '$%s$ & %.2f & %.2f & %s \\\\\n', ...
        FtestResults.Variable(i), ...
        FtestResults.F_stat(i), ...
        FtestResults.p_value(i), ...
        FtestResults.Decision(i));
end

% LaTeX table footer
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{F-test results for $H_0: \\beta_2 + \\gamma_2 = 0$ in each regression.}\n');
fprintf(fid, '\\label{tab:Ftest_results}\n');

% Close the file
fclose(fid);

disp('LaTeX table exported!');



