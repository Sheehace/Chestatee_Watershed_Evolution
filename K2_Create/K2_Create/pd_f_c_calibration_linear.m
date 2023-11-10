clear all
close all
clc



%% Variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_type = 'png'
file_res = 600

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create Output_pd_f_c_calibration_linear folder
mkdir('Output_pd_f_c_calibration_linear');





%%

% Create Output folder
mkdir('Output');

% Import DEM
script_directory = pwd;
DEM_path = strrep(script_directory, 'K2_Create', 'DEM/Chestatee.asc');
DEM = GRIDobj(DEM_path);

% Import Coords
Coords = readmatrix('Coords.csv');

% Plot
fig = figure(1);
hold on;
imageschs(DEM, [], 'colormap', 'turbo');
plot(Coords(:, 7), Coords(:, 6), 'ko', 'MarkerFaceColor', 'r');
exportgraphics(fig, ['Output/MACA_Input_Coordinate_Map.', file_type], 'Resolution', file_res)
close fig 1
%
fig = figure(1);
hold on;
imageschs(DEM, [], 'colormap', 'turbo');
plot(Coords(:, 8), Coords(:, 9), 'ks', 'MarkerFaceColor', 'y');
exportgraphics(fig, ['Output/MACA_Cell_Center_Map.', file_type], 'Resolution', file_res)
close fig 1
%
fig = figure(1);
hold on;
imageschs(DEM, [], 'colormap', 'turbo');
plot(Coords(:, 7), Coords(:, 6), 'ko', 'MarkerFaceColor', 'r');
plot(Coords(:, 8), Coords(:, 9), 'ks', 'MarkerFaceColor', 'y');
exportgraphics(fig, ['Output/MACA_Combined_Map.', file_type], 'Resolution', file_res)
close fig 1




%% Hist

% Create "List" of csv files
% folder = 'C:/Users/sheehacz/Dropbox/BC_Landlab/SPACE_2022/New_Models/Calib/Climate/Models/Linear/Hist';
folder = 'MACA_Data/Linear/Hist';
contents = dir(folder);
names = {contents.name};
List = {};
i = 1;
for c = 1 : size(names, 2)
    if contains(names{c}, '.csv') == 1 
        List{i} = names{c};
        i = i + 1;       
    end
end
%
clear c;
clear contents;
clear i;
clear names;

% Initiate "C"
Pd = (1 : size(List, 2)) * NaN;
F = (1 : size(List, 2)) * NaN;
C = (1 : size(List, 2)) * NaN;

% Create Output_pd_f_c_calibration_linear folder
mkdir('Output_pd_f_c_calibration_linear/Hist');

% Add csv files to "p"
for i = 1 : size(List, 2)
    
    % Create Output_pd_f_c_calibration_linear folder
    mkdir(['Output_pd_f_c_calibration_linear/Hist/' num2str(i)])

    % Import file
    p = readmatrix([folder '/' num2str(i) '.csv']);
    
    % Calculate record length (days)
    record_length = size(p, 1);
    
    % Delete time column
    p(:, 1) = [];
    
    % Remove nans and zeros
    p(isnan(p)) = [];       % Remove nans
    p(p == 0) = [];         % Remove zeros

    % Calculate "pd" (method of doing this partially inside loop above and in 
    % this line is for consistancy with Python script that calculates pd and F).
    pd = mean(p);
    Pd(i) = pd;
    
    % Calculate f
    f = size(p, 1) / record_length;
    F(i) = f;
    
    % Find unique values of p
    p_unique = unique(p);

    % Calculate number of events larger than each p_unique
    count = p_unique * nan;
    for j = 1 : size(p_unique, 1)
        count(j) = sum(p > p_unique(j));  
    end

    % Calculate exceedence frequency
    Pr = count / size(p, 1);

    % Plot Pr
    fig1 = figure(1);
    semilogy(p_unique, Pr, '.');
    ax = gca;
    color = ax.ColorOrder(ax.ColorOrderIndex - 1, :);
    
    % Calculate double ln transformed data
    y = log(Pr);
    y = y * -1;           
    y = log(y);
    %
    x = log(p_unique);
    y(end) = [];
    x(end) = [];
    
%     % Find datapoints over 95th percentile and Fit data to find "c" 
%     perc95 = prctile(Pr, 5);
%     include = Pr >= perc95;
%     size_include = size(Pr(include), 1);
%     pf = polyfit(x(include), y(include), 1);
%     pfy = polyval(pf, x(include));
%     c = pf(1);
    pf = polyfit(x(x>0), y(x>0), 1);
    c = pf(1);
    
    % Store "c"
    C(i) = c;
    
    % Plot log-transformed data and export
    fig2 = figure(2);
    hold on;
    plot(x, y, '.');
    plot(x, (x*c) + pf(2));
    xlabel('ln(p)');
    ylabel('ln( ln(Pr) )');
    legend('data', 'fit');
    txt = ['c = ', num2str(c)];
    text(4, -6, txt);
    exportgraphics(fig2, ['Output_pd_f_c_calibration_linear/Hist/' num2str(i) '/Log_Transformed.' file_type], 'Resolution', file_res)
    close fig 2;
    
    % Add extrapolation to Figure 1 and export
    figure(1);
    hold on;
    lamda = pd / (gamma(1 + (1 / c)));
    fit = exp( -(p_unique / lamda) .^ c );
    plot(p_unique, fit, 'Color', color);
    xlabel('p');
    ylabel('Pr');
    legend('data', 'fit');
    txt = ['c = ', num2str(c)];
    text(20, 10^-2, txt);
    exportgraphics(fig1, ['Output_pd_f_c_calibration_linear/Hist/' num2str(i) '/p_vs_ExhedanceFrequancy.' file_type], 'Resolution', file_res)
    close fig 1;
    
end

% Store data
C_Hist_mean = mean(C);
C_Hist_std = std(C);
Pd_Hist_mean = mean(Pd);
Pd_Hist_std = std(Pd);
F_Hist_mean = mean(F);
F_Hist_std = std(F);





%% RCP 4.5

% Create "List" of csv files
% folder = 'C:/Users/sheehacz/Dropbox/BC_Landlab/SPACE_2022/New_Models/Calib/Climate/Models/Linear/RCP45';
folder = 'MACA_Data/Linear/RCP45';
contents = dir(folder);
names = {contents.name};
List = {};
i = 1;
for c = 1 : size(names, 2)
    if contains(names{c}, '.csv') == 1 
        List{i} = names{c};
        i = i + 1;       
    end
end
%
clear c;
clear contents;
clear i;
clear names;

% Initiate "C"
Pd = (1 : size(List, 2)) * NaN;
F = (1 : size(List, 2)) * NaN;
C = (1 : size(List, 2)) * NaN;

% Create Output_pd_f_c_calibration_linear folder
mkdir('Output_pd_f_c_calibration_linear/RCP45');

% Add csv files to "p"
for i = 1 : size(List, 2)
    
    % Create Output_pd_f_c_calibration_linear folder
    mkdir(['Output_pd_f_c_calibration_linear/RCP45/' num2str(i)])

    % Import file
    p = readmatrix([folder '/' num2str(i) '.csv']);
    
    % Calculate record length (days)
    record_length = size(p, 1);
    
    % Delete time column
    p(:, 1) = [];
    
    % Remove nans and zeros
    p(isnan(p)) = [];       % Remove nans
    p(p == 0) = [];         % Remove zeros

    % Calculate "pd" (method of doing this partially inside loop above and in 
    % this line is for consistancy with Python script that calculates pd and F).
    pd = mean(p);
    Pd(i) = pd;
    
    % Calculate f
    f = size(p, 1) / record_length;
    F(i) = f;
    
    % Find unique values of p
    p_unique = unique(p);

    % Calculate number of events larger than each p_unique
    count = p_unique * nan;
    for j = 1 : size(p_unique, 1)
        count(j) = sum(p > p_unique(j));  
    end

    % Calculate exceedence frequency
    Pr = count / size(p, 1);

    % Plot Pr
    fig1 = figure(1);
    semilogy(p_unique, Pr, '.');
    ax = gca;
    color = ax.ColorOrder(ax.ColorOrderIndex - 1, :);
    
    % Calculate double ln transformed data
    y = log(Pr);
    y = y * -1;           
    y = log(y);
    %
    x = log(p_unique);
    y(end) = [];
    x(end) = [];
    
%     % Find datapoints over 95th percentile and Fit data to find "c" 
%     perc95 = prctile(Pr, 5);
%     include = Pr >= perc95;
%     size_include = size(Pr(include), 1);
%     pf = polyfit(x(include), y(include), 1);
%     pfy = polyval(pf, x(include));
%     c = pf(1);
    pf = polyfit(x(x>0), y(x>0), 1);
    c = pf(1);
    
    % Store "c"
    C(i) = c;
    
    % Plot log-transformed data and export
    fig2 = figure(2);
    hold on;
    plot(x, y, '.');
    plot(x, (x*c) + pf(2));
    xlabel('ln(p)');
    ylabel('ln( ln(Pr) )');
    legend('data', 'fit');
    txt = ['c = ', num2str(c)];
    text(4, -6, txt);
    exportgraphics(fig2, ['Output_pd_f_c_calibration_linear/RCP45/' num2str(i) '/Log_Transformed.' file_type], 'Resolution', file_res)
    close fig 2;
    
    % Add extrapolation to Figure 1 and export
    figure(1);
    hold on;
    lamda = pd / (gamma(1 + (1 / c)));
    fit = exp( -(p_unique / lamda) .^ c );
    plot(p_unique, fit, 'Color', color);
    xlabel('p');
    ylabel('Pr');
    legend('data', 'fit');
    txt = ['c = ', num2str(c)];
    text(20, 10^-2, txt);
    exportgraphics(fig1, ['Output_pd_f_c_calibration_linear/RCP45/' num2str(i) '/p_vs_ExhedanceFrequancy.' file_type], 'Resolution', file_res)
    close fig 1;
    
end

% Store data
C_RCP45_mean = mean(C);
C_RCP45_std = std(C);
Pd_RCP45_mean = mean(Pd);
Pd_RCP45_std = std(Pd);
F_RCP45_mean = mean(F);
F_RCP45_std = std(F);





%% RCP 8.5

% Create "List" of csv files
% folder = 'C:/Users/sheehacz/Dropbox/BC_Landlab/SPACE_2022/New_Models/Calib/Climate/Models/Linear/RCP85';
folder = 'MACA_Data/Linear/RCP85';
contents = dir(folder);
names = {contents.name};
List = {};
i = 1;
for c = 1 : size(names, 2)
    if contains(names{c}, '.csv') == 1 
        List{i} = names{c};
        i = i + 1;       
    end
end
%
clear c;
clear contents;
clear i;
clear names;

% Initiate "C"
Pd = (1 : size(List, 2)) * NaN;
F = (1 : size(List, 2)) * NaN;
C = (1 : size(List, 2)) * NaN;

% Create Output_pd_f_c_calibration_linear folder
mkdir('Output_pd_f_c_calibration_linear/RCP85');

% Add csv files to "p"
for i = 1 : size(List, 2)
    
    % Create Output_pd_f_c_calibration_linear folder
    mkdir(['Output_pd_f_c_calibration_linear/RCP85/' num2str(i)])

    % Import file
    p = readmatrix([folder '/' num2str(i) '.csv']);
    
    % Calculate record length (days)
    record_length = size(p, 1);
    
    % Delete time column
    p(:, 1) = [];
    
    % Remove nans and zeros
    p(isnan(p)) = [];       % Remove nans
    p(p == 0) = [];         % Remove zeros

    % Calculate "pd" (method of doing this partially inside loop above and in 
    % this line is for consistancy with Python script that calculates pd and F).
    pd = mean(p);
    Pd(i) = pd;
    
    % Calculate f
    f = size(p, 1) / record_length;
    F(i) = f;
    
    % Find unique values of p
    p_unique = unique(p);

    % Calculate number of events larger than each p_unique
    count = p_unique * nan;
    for j = 1 : size(p_unique, 1)
        count(j) = sum(p > p_unique(j));  
    end

    % Calculate exceedence frequency
    Pr = count / size(p, 1);

    % Plot Pr
    fig1 = figure(1);
    semilogy(p_unique, Pr, '.');
    ax = gca;
    color = ax.ColorOrder(ax.ColorOrderIndex - 1, :);
    
    % Calculate double ln transformed data
    y = log(Pr);
    y = y * -1;           
    y = log(y);
    %
    x = log(p_unique);
    y(end) = [];
    x(end) = [];
    
%     % Find datapoints over 95th percentile and Fit data to find "c" 
%     perc95 = prctile(Pr, 5);
%     include = Pr >= perc95;
%     size_include = size(Pr(include), 1);
%     pf = polyfit(x(include), y(include), 1);
%     pfy = polyval(pf, x(include));
%     c = pf(1);
    pf = polyfit(x(x>0), y(x>0), 1);
    c = pf(1);
    
    % Store "c"
    C(i) = c;
    
    % Plot log-transformed data and export
    fig2 = figure(2);
    hold on;
    plot(x, y, '.');
    plot(x, (x*c) + pf(2));
    xlabel('ln(p)');
    ylabel('ln( ln(Pr) )');
    legend('data', 'fit');
    txt = ['c = ', num2str(c)];
    text(4, -6, txt);
    exportgraphics(fig2, ['Output_pd_f_c_calibration_linear/RCP85/' num2str(i) '/Log_Transformed.' file_type], 'Resolution', file_res)
    close fig 2;
    
    % Add extrapolation to Figure 1 and export
    figure(1);
    hold on;
    lamda = pd / (gamma(1 + (1 / c)));
    fit = exp( -(p_unique / lamda) .^ c );
    plot(p_unique, fit, 'Color', color);
    xlabel('p');
    ylabel('Pr');
    legend('data', 'fit');
    txt = ['c = ', num2str(c)];
    text(20, 10^-2, txt);
    exportgraphics(fig1, ['Output_pd_f_c_calibration_linear/RCP85/' num2str(i) '/p_vs_ExhedanceFrequancy.' file_type], 'Resolution', file_res)
    close fig 1;
    
end

% Store data
C_RCP85_mean = mean(C);
C_RCP85_std = std(C);
Pd_RCP85_mean = mean(Pd);
Pd_RCP85_std = std(Pd);
F_RCP85_mean = mean(F);
F_RCP85_std = std(F);







%% Following sections brought over from "K_Ratio.m" script 

%% Set up

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dp = 0.001
pmax = 1000
Im = 15
m = 0.503

file_type = 'png'
file_res = 600

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Calculate RCP45_Ratio

% Assign variables
pd0 = Pd_Hist_mean;
pd1 = Pd_RCP45_mean;
F0 = F_Hist_mean;
F1 = F_RCP45_mean;
c0 = C_Hist_mean;
c1 = C_RCP45_mean;

% Set p array
p = 0.01 : dp : pmax;

% Calculate lamdas
lamda0 = pd0 / gamma(1 + (1 / c0));
lamda1 = pd1 / gamma(1 + (1 / c1));

% Set index
index = 1;

% Integrate
for j = Im : dp : pmax

    % Calc numerator and denominator
    Den(index) = ((j - Im)^m) * ( (c0 / lamda0) .* ( (j / lamda0) .^ (c0 - 1) ) .* ( exp( -(j / lamda0) .^ c0 ) ) ) * dp;
    Num(index) = ((j - Im)^m) * ( (c1 / lamda1) .* ( (j / lamda1) .^ (c1 - 1) ) .* ( exp( -(j / lamda1) .^ c1 ) ) ) * dp;

    % Advance index
    index = index + 1;

end

% Fractionalize
Den = Den * F0;
Num = Num * F1;

%Integrate
Num_sum = nansum(Num);
Den_sum = nansum(Den);

% Divide
RCP45_Ratio = Num_sum / Den_sum;





%% Calculate RCP85_Ratio

% Assign variables
pd0 = Pd_Hist_mean;
pd1 = Pd_RCP85_mean;
F0 = F_Hist_mean;
F1 = F_RCP85_mean;
c0 = C_Hist_mean;
c1 = C_RCP85_mean;

% Set p array
p = 0.01 : dp : pmax;

% Calculate lamdas
lamda0 = pd0 / gamma(1 + (1 / c0));
lamda1 = pd1 / gamma(1 + (1 / c1));

% Set index
index = 1;

% Integrate
for j = Im : dp : pmax

    % Calc numerator and denominator
    Den(index) = ((j - Im)^m) * ( (c0 / lamda0) .* ( (j / lamda0) .^ (c0 - 1) ) .* ( exp( -(j / lamda0) .^ c0 ) ) ) * dp;
    Num(index) = ((j - Im)^m) * ( (c1 / lamda1) .* ( (j / lamda1) .^ (c1 - 1) ) .* ( exp( -(j / lamda1) .^ c1 ) ) ) * dp;

    % Advance index
    index = index + 1;

end

% Fractionalize
Den = Den * F0;
Num = Num * F1;

%Integrate
Num_sum = nansum(Num);
Den_sum = nansum(Den);

% Divide
RCP85_Ratio = Num_sum / Den_sum;





%% Write csv

% Table
T = table(RCP45_Ratio, RCP85_Ratio, C_Hist_mean, C_Hist_std, C_RCP45_mean, C_RCP45_std, C_RCP85_mean, C_RCP85_std, Pd_Hist_mean, Pd_Hist_std, Pd_RCP45_mean, Pd_RCP45_std, Pd_RCP85_mean, Pd_RCP85_std, F_Hist_mean, F_Hist_std, F_RCP45_mean, F_RCP45_std, F_RCP85_mean, F_RCP85_std);
writetable(T, 'Output_pd_f_c_calibration_linear/Data.csv');




