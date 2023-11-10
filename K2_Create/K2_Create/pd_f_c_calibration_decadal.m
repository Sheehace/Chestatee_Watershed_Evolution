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

% Create output folder
mkdir('Output_pd_f_c_calibration_decadal');





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





%% Hist (copied from pd_f_c_calibration_linear.m)

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
mkdir('Output_pd_f_c_calibration_decadal/Hist');

% Add csv files to "p"
for i = 1 : size(List, 2)
    
    % Create Output_pd_f_c_calibration_linear folder
    mkdir(['Output_pd_f_c_calibration_decadal/Hist/' num2str(i)])

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
    exportgraphics(fig2, ['Output_pd_f_c_calibration_decadal/Hist/' num2str(i) '/Log_Transformed.' file_type], 'Resolution', file_res)
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
    exportgraphics(fig1, ['Output_pd_f_c_calibration_decadal/Hist/' num2str(i) '/p_vs_ExhedanceFrequancy.' file_type], 'Resolution', file_res)
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
% folder = 'C:/Users/sheehacz/Dropbox/BC_Landlab/SPACE_2022/New_Models/Calib/Climate/Models/Decadal/RCP45';
folder = 'MACA_Data/Decadal/Combined/RCP45';
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

% Create output folder
mkdir('Output_pd_f_c_calibration_decadal/RCP45');

% Add csv files to "p"
for i = 1 : size(List, 2)
    
    % Create output folder
    mkdir(['Output_pd_f_c_calibration_decadal/RCP45/' num2str(i)])

    % Import file
    %p = readmatrix([folder '/' num2str(i) '.csv']);
    p = readtable([folder '/' num2str(i) '.csv']);
    p.Properties.VariableNames = ["Date", "Prcp"];
    
    % Extract years
    Years = year(p.Date);
    
    % Initiate starting year index
    year_lower = 2000;
    year_upper = 2010;
    
    % Set indeces for loop below
    toggle = 1;
    row = 1;
    
    % Loop through decades
    while toggle == 1
        
        % Create directory for current decade
        mkdir(['Output_pd_f_c_calibration_decadal/RCP45/' num2str(i) '/' num2str(year_lower)]);
        
        % Identify days within current decade
        use = Years >= year_lower & Years < year_upper;
        
        % Create subset of "p" only containing relevant decade
        subset = p(use, :);
        
        % Calculate record length (days)
        record_length = size(subset, 1);
        
        % Copy rainfall to separate array
        %subset(:, 1) = [];
        prcp = subset.Prcp;
        
        % Remove nans and zeros
        prcp(isnan(prcp)) = [];       % Remove nans
        prcp(prcp == 0) = [];         % Remove zeros
        
        % Calculate "pd" (method of doing this partially inside loop above and in 
        % this line is for consistancy with Python script that calculates pd and F).
        pd = mean(prcp);
        Pd(row, i) = pd;
        
        % Calculate f
        f = size(prcp, 1) / record_length;
        F(row, i) = f;
        
        % Find unique values of p
        prcp_unique = unique(prcp);
        
        % Calculate number of events larger than each p_unique
        count = prcp_unique * nan;
        for j = 1 : size(prcp_unique, 1)
            count(j) = sum(subset.Prcp > prcp_unique(j));  
        end
        
        % Calculate exceedence frequency
        Pr = count / size(prcp, 1);
        
        % Plot Pr
        fig1 = figure(1);
        semilogy(prcp_unique, Pr, '.');
        ax = gca;
        %color = ax.ColorOrder(ax.ColorOrderIndex - 1, :);
    
        % Calculate double ln transformed data
        y = log(Pr);
        y = y * -1;           
        y = log(y);
        %
        x = log(prcp_unique);
        y(end) = [];
        x(end) = [];
        %
        
%         % Find datapoints over 95th percentile and Fit data to find "c" 
%         perc95 = prctile(Pr, 5);
%         include = Pr >= perc95;
%         size_include = size(Pr(include), 1);
%         pf = polyfit(x(include), y(include), 1);
%         pfy = polyval(pf, x(include));
%         c = pf(1);        
        pf = polyfit(x(x>0), y(x>0), 1);
        c = pf(1);
        
        % Store "c"
        C(row, i) = c;
        
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
        exportgraphics(fig2, ['Output_pd_f_c_calibration_decadal/RCP45/' num2str(i) '/' num2str(year_lower) '/Log_Transformed.' file_type], 'Resolution', file_res)
        close fig 2;
        
        % Add extrapolation to Figure 1 and export
        figure(1);
        hold on;
        lamda = pd / (gamma(1 + (1 / c)));
        fit = exp( -(prcp_unique / lamda) .^ c );
        %plot(prcp_unique, fit, 'Color', color);
        plot(prcp_unique, fit);
        xlabel('p');
        ylabel('Pr');
        legend('data', 'fit');
        txt = ['c = ', num2str(c)];
        text(20, 10^-2, txt);
        exportgraphics(fig1, ['Output_pd_f_c_calibration_decadal/RCP45/' num2str(i) '/' num2str(year_lower) '/p_vs_ExhedanceFrequancy.' file_type], 'Resolution', file_res)
        close fig 1;
        
        % Advance indeces
        year_lower = year_lower + 10;
        year_upper = year_upper + 10;
        row = row + 1;
        
        if year_upper == 2110
            toggle = 0;
        end
  
    end
    
end

% Store data
C_RCP45_mean = NaN;
C_RCP45_std = NaN;
Pd_RCP45_mean = NaN;
Pd_RCP45_std = NaN;
F_RCP45_mean = NaN;
F_RCP45_std = NaN;
%
for r = 1 : size(Pd, 1)
    
    C_RCP45_mean(r, 1) = mean(C(r, :));
    C_RCP45_std(r,1) = std(C(r, :));
    Pd_RCP45_mean(r, 1) = mean(Pd(r, :));
    Pd_RCP45_std(r, 1) = std(Pd(r, :));
    F_RCP45_mean(r, 1) = mean(F(r, :));
    F_RCP45_std(r, 1) = std(F(r, :));

end

% Plot Pd
fig = figure();
x = 2000 : 10 : 2090;
errorbar(x, Pd_RCP45_mean(:, 1), Pd_RCP45_std(:, 1), '.-');
xlabel('Decade');
ylabel('Pd');
exportgraphics(fig, ['Output_pd_f_c_calibration_decadal/RCP45/Pd.'  file_type], 'Resolution', file_res);
close fig 1;

% Plot F
fig = figure();
x = 2000 : 10 : 2090;
errorbar(x, F_RCP45_mean(:, 1), F_RCP45_std(:, 1), '.-');
xlabel('Decade');
ylabel('F');
exportgraphics(fig, ['Output_pd_f_c_calibration_decadal/RCP45/F.'  file_type], 'Resolution', file_res);
close fig 1;

% Plot F
fig = figure();
x = 2000 : 10 : 2090;
errorbar(x, C_RCP45_mean(:, 1), C_RCP45_std(:, 1), '.-');
xlabel('Decade');
ylabel('C');
exportgraphics(fig, ['Output_pd_f_c_calibration_decadal/RCP45/C.'  file_type], 'Resolution', file_res);
close fig 1;





%% RCP 8.5

% Create "List" of csv files
% folder = 'C:/Users/sheehacz/Dropbox/BC_Landlab/SPACE_2022/New_Models/Calib/Climate/Models/Decadal/RCP85';
folder = 'MACA_Data/Decadal/Combined/RCP85';
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

% Create output folder
mkdir('Output_pd_f_c_calibration_decadal/RCP85');

% Add csv files to "p"
for i = 1 : size(List, 2)
    
    % Create output folder
    mkdir(['Output_pd_f_c_calibration_decadal/RCP85/' num2str(i)])

    % Import file
    %p = readmatrix([folder '/' num2str(i) '.csv']);
    p = readtable([folder '/' num2str(i) '.csv']);
    p.Properties.VariableNames = ["Date", "Prcp"];
    
    % Extract years
    Years = year(p.Date);
    
    % Initiate starting year index
    year_lower = 2000;
    year_upper = 2010;
    
    % Set indeces for loop below
    toggle = 1;
    row = 1;
    
    % Loop through decades
    while toggle == 1
        
        % Create directory for current decade
        mkdir(['Output_pd_f_c_calibration_decadal/RCP85/' num2str(i) '/' num2str(year_lower)]);
        
        % Identify days within current decade
        use = Years >= year_lower & Years < year_upper;
        
        % Create subset of "p" only containing relevant decade
        subset = p(use, :);
        
        % Calculate record length (days)
        record_length = size(subset, 1);
        
        % Copy rainfall to separate array
        %subset(:, 1) = [];
        prcp = subset.Prcp;
        
        % Remove nans and zeros
        prcp(isnan(prcp)) = [];       % Remove nans
        prcp(prcp == 0) = [];         % Remove zeros
        
        % Calculate "pd" (method of doing this partially inside loop above and in 
        % this line is for consistancy with Python script that calculates pd and F).
        pd = mean(prcp);
        Pd(row, i) = pd;
        
        % Calculate f
        f = size(prcp, 1) / record_length;
        F(row, i) = f;
        
        % Find unique values of p
        prcp_unique = unique(prcp);
        
        % Calculate number of events larger than each p_unique
        count = prcp_unique * nan;
        for j = 1 : size(prcp_unique, 1)
            count(j) = sum(subset.Prcp > prcp_unique(j));  
        end
        
        % Calculate exceedence frequency
        Pr = count / size(prcp, 1);
        
        % Plot Pr
        fig1 = figure(1);
        semilogy(prcp_unique, Pr, '.');
        ax = gca;
        %color = ax.ColorOrder(ax.ColorOrderIndex - 1, :);
    
        % Calculate double ln transformed data
        y = log(Pr);
        y = y * -1;           
        y = log(y);
        %
        x = log(prcp_unique);
        y(end) = [];
        x(end) = [];
        %
        
%         % Find datapoints over 95th percentile and Fit data to find "c" 
%         perc95 = prctile(Pr, 5);
%         include = Pr >= perc95;
%         size_include = size(Pr(include), 1);
%         pf = polyfit(x(include), y(include), 1);
%         pfy = polyval(pf, x(include));
%         c = pf(1);        
        pf = polyfit(x(x>0), y(x>0), 1);
        c = pf(1);
        
        % Store "c"
        C(row, i) = c;
        
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
        exportgraphics(fig2, ['Output_pd_f_c_calibration_decadal/RCP85/' num2str(i) '/' num2str(year_lower) '/Log_Transformed.' file_type], 'Resolution', file_res)
        close fig 2;
        
        % Add extrapolation to Figure 1 and export
        figure(1);
        hold on;
        lamda = pd / (gamma(1 + (1 / c)));
        fit = exp( -(prcp_unique / lamda) .^ c );
        %plot(prcp_unique, fit, 'Color', color);
        plot(prcp_unique, fit);
        xlabel('p');
        ylabel('Pr');
        legend('data', 'fit');
        txt = ['c = ', num2str(c)];
        text(20, 10^-2, txt);
        exportgraphics(fig1, ['Output_pd_f_c_calibration_decadal/RCP85/' num2str(i) '/' num2str(year_lower) '/p_vs_ExhedanceFrequancy.' file_type], 'Resolution', file_res)
        close fig 1;
        
        % Advance indeces
        year_lower = year_lower + 10;
        year_upper = year_upper + 10;
        row = row + 1;
        
        if year_upper == 2110
            toggle = 0;
        end
  
    end
    
end

% Store data
C_RCP85_mean = NaN;
C_RCP85_std = NaN;
Pd_RCP85_mean = NaN;
Pd_RCP85_std = NaN;
F_RCP85_mean = NaN;
F_RCP85_std = NaN;
%
for r = 1 : size(Pd, 1)
    
    C_RCP85_mean(r, 1) = mean(C(r, :));
    C_RCP85_std(r,1) = std(C(r, :));
    Pd_RCP85_mean(r, 1) = mean(Pd(r, :));
    Pd_RCP85_std(r, 1) = std(Pd(r, :));
    F_RCP85_mean(r, 1) = mean(F(r, :));
    F_RCP85_std(r, 1) = std(F(r, :));

end

% Plot Pd
fig = figure();
x = 2000 : 10 : 2090;
errorbar(x, Pd_RCP85_mean(:, 1), Pd_RCP85_std(:, 1), '.-');
xlabel('Decade');
ylabel('Pd');
exportgraphics(fig, ['Output_pd_f_c_calibration_decadal/RCP85/Pd.'  file_type], 'Resolution', file_res);
close fig 1;

% Plot F
fig = figure();
x = 2000 : 10 : 2090;
errorbar(x, F_RCP85_mean(:, 1), F_RCP85_std(:, 1), '.-');
xlabel('Decade');
ylabel('F');
exportgraphics(fig, ['Output_pd_f_c_calibration_decadal/RCP85/F.'  file_type], 'Resolution', file_res);
close fig 1;

% Plot F
fig = figure();
x = 2000 : 10 : 2090;
errorbar(x, C_RCP85_mean(:, 1), C_RCP85_std(:, 1), '.-');
xlabel('Decade');
ylabel('C');
exportgraphics(fig, ['Output_pd_f_c_calibration_decadal/RCP85/C.'  file_type], 'Resolution', file_res);
close fig 1;


% FINISHED HERE 2022 10 17 Before Washington DC





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

RCP45_Ratio = NaN;

for r = 1 : size(Pd_RCP45_mean, 1)

    % Assign variables
    %pd0 = Pd_RCP45_mean(1);
    %pd0 = [9.22296461023816];   % Value calculated in linear code for historical 1950 - 2000
    pd0 = Pd_Hist_mean;
    pd1 = Pd_RCP45_mean(r);
    %F0 = F_RCP45_mean(1);
    %F0 = [0.477209284166553];   % Value calculated in linear code for historical 1950 - 2000
    F0 = F_Hist_mean;
    F1 = F_RCP45_mean(r);
    %c0 = C_RCP45_mean(1);
    %c0 = [0.716885884283468];   % Value calculated in linear code for historical 1950 - 2000
    c0 = C_Hist_mean;
    c1 = C_RCP45_mean(r);

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
    RCP45_Ratio(r, 1) = Num_sum / Den_sum;


end

% Plot K ratio
x = 2000 : 10 : 2090;
fig = figure();
hold on;
plot(x, RCP45_Ratio, 'o-');
xlabel('Decade');
ylabel('K Ratio');
title('RCP 4.5');
exportgraphics(fig, ['Output_pd_f_c_calibration_decadal/RCP45.'  file_type], 'Resolution', file_res);
close fig 1;





%% Calculate RCP85_Ratio

RCP85_Ratio = NaN;

for r = 1 : size(Pd_RCP85_mean, 1)

    % Assign variables
    %pd0 = Pd_RCP85_mean(1);
    %pd0 = [9.22296461023816];   % Value calculated in linear code for historical 1950 - 2000
    pd0 = Pd_Hist_mean;
    pd1 = Pd_RCP85_mean(r);
    %F0 = F_RCP45_mean(1);
    %F0 = [0.477209284166553];   % Value calculated in linear code for historical 1950 - 2000
    F0 = F_Hist_mean;
    F1 = F_RCP85_mean(r);
    %c0 = C_RCP45_mean(1);
    %c0 = [0.716885884283468];   % Value calculated in linear code for historical 1950 - 2000
    c0 = C_Hist_mean;
    c1 = C_RCP85_mean(r);

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

    % Integrate
    Num_sum = nansum(Num);
    Den_sum = nansum(Den);

    % Divide
    RCP85_Ratio(r, 1) = Num_sum / Den_sum;


end

% Plot K ratio
x = 2000 : 10 : 2090;
fig = figure();
hold on;
plot(x, RCP85_Ratio, 'o-');
xlabel('Decade');
ylabel('K Ratio');
title('RCP 8.5');
exportgraphics(fig, ['Output_pd_f_c_calibration_decadal/RCP85.'  file_type], 'Resolution', file_res);
close fig 1;





%% Plot combined and write csv

% Plot combined K ratio
x = 2000 : 10 : 2090;
fig = figure();
hold on;
plot(x, RCP45_Ratio, 'o-');
plot(x, RCP85_Ratio, 'o-');
xlabel('Decade');
ylabel('K Ratio');
title('RCP 8.5');
legend('RCP 4.5', 'RCP 8.5')
exportgraphics(fig, ['Output_pd_f_c_calibration_decadal/K_Ratio_Combined.'  file_type], 'Resolution', file_res);
close fig 1;

% Table
T45 = table(RCP45_Ratio);
T85 = table(RCP85_Ratio);
writetable(T45, 'Output_pd_f_c_calibration_decadal/RCP45_Data.csv');
writetable(T85, 'Output_pd_f_c_calibration_decadal/RCP85_Data.csv');





%% Calculate RCP45_Ratio lower and upper limit

% Lower Limit

RCP45_Ratio_low = NaN;

for r = 1 : size(Pd_RCP45_mean, 1)

    % Assign variables
    %pd0 = Pd_RCP45_mean(1);
    %pd0 = [9.22296461023816];   % Value calculated in linear code for historical 1950 - 2000
    pd0 = Pd_Hist_mean;
    pd1 = Pd_RCP45_mean(r) - Pd_RCP45_std(r);
    %F0 = F_RCP45_mean(1);
    %F0 = [0.477209284166553];   % Value calculated in linear code for historical 1950 - 2000
    F0 = F_Hist_mean;
    F1 = F_RCP45_mean(r) - F_RCP45_std(r);
    %c0 = C_RCP45_mean(1);
    %c0 = [0.716885884283468];   % Value calculated in linear code for historical 1950 - 2000
    c0 = C_Hist_mean;
    c1 = C_RCP45_mean(r) + C_RCP45_std(r);

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
    RCP45_Ratio_low(r, 1) = Num_sum / Den_sum;


end

% Upper Limit

RCP45_Ratio_high = NaN;

for r = 1 : size(Pd_RCP45_mean, 1)

    % Assign variables
    %pd0 = Pd_RCP45_mean(1);
    %pd0 = [9.22296461023816];   % Value calculated in linear code for historical 1950 - 2000
    pd0 = Pd_Hist_mean;
    pd1 = Pd_RCP45_mean(r) + Pd_RCP45_std(r);
    %F0 = F_RCP45_mean(1);
    %F0 = [0.477209284166553];   % Value calculated in linear code for historical 1950 - 2000
    F0 = F_Hist_mean;
    F1 = F_RCP45_mean(r) + F_RCP45_std(r);
    %c0 = C_RCP45_mean(1);
    %c0 = [0.716885884283468];   % Value calculated in linear code for historical 1950 - 2000
    c0 = C_Hist_mean;
    c1 = C_RCP45_mean(r) - C_RCP45_std(r);

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
    RCP45_Ratio_high(r, 1) = Num_sum / Den_sum;


end

% Plot K ratio w/ stds
x = 2000 : 10 : 2090;
neg = RCP45_Ratio - RCP45_Ratio_low;
pos = RCP45_Ratio_high - RCP45_Ratio;
fig = figure();
hold on;
errorbar(x, RCP45_Ratio, neg, pos, 'o-');
xlabel('Decade');
ylabel('K Ratio');
title('RCP 4.5');
exportgraphics(fig, ['Output_pd_f_c_calibration_decadal/K_Ratio_45_errorbar.'  file_type], 'Resolution', file_res);
close fig 1;





%% Calculate RCP85_Ratio lower and upper limit

% Lower Limit

RCP85_Ratio_low = NaN;

for r = 1 : size(Pd_RCP85_mean, 1)

    % Assign variables
    %pd0 = Pd_RCP85_mean(1);
    %pd0 = [9.22296461023816];   % Value calculated in linear code for historical 1950 - 2000
    pd0 = Pd_Hist_mean;
    pd1 = Pd_RCP85_mean(r) - Pd_RCP85_std(r);
    %F0 = F_RCP85_mean(1);
    %F0 = [0.477209284166553];   % Value calculated in linear code for historical 1950 - 2000
    F0 = F_Hist_mean;
    F1 = F_RCP85_mean(r) - F_RCP85_std(r);
    %c0 = C_RCP85_mean(1);
    %c0 = [0.716885884283468];   % Value calculated in linear code for historical 1950 - 2000
    c0 = C_Hist_mean;
    c1 = C_RCP85_mean(r) + C_RCP85_std(r);

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
    RCP85_Ratio_low(r, 1) = Num_sum / Den_sum;


end

% Upper Limit

RCP85_Ratio_high = NaN;

for r = 1 : size(Pd_RCP85_mean, 1)

    % Assign variables
    %pd0 = Pd_RCP85_mean(1);
    %pd0 = [9.22296461023816];   % Value calculated in linear code for historical 1950 - 2000
    pd0 = Pd_Hist_mean;
    pd1 = Pd_RCP85_mean(r) + Pd_RCP85_std(r);
    %F0 = F_RCP85_mean(1);
    %F0 = [0.477209284166553];   % Value calculated in linear code for historical 1950 - 2000
    F0 = F_Hist_mean;
    F1 = F_RCP85_mean(r) + F_RCP85_std(r);
    %c0 = C_RCP85_mean(1);
    %c0 = [0.716885884283468];   % Value calculated in linear code for historical 1950 - 2000
    c0 = C_Hist_mean;
    c1 = C_RCP85_mean(r) - C_RCP85_std(r);

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
    RCP85_Ratio_high(r, 1) = Num_sum / Den_sum;


end

% Plot K ratio w/ stds
x = 2000 : 10 : 2090;
neg = RCP85_Ratio - RCP85_Ratio_low;
pos = RCP85_Ratio_high - RCP85_Ratio;
fig = figure();
hold on;
errorbar(x, RCP85_Ratio, neg, pos, 'o-');
xlabel('Decade');
ylabel('K Ratio');
title('RCP 8.5');
exportgraphics(fig, ['Output_pd_f_c_calibration_decadal/K_Ratio_85_errorbar.'  file_type], 'Resolution', file_res);
close fig 1;





%% Plot combined errorbars

% Pd_Hist_mean = [9.22296461023816];
% F_Hist_mean = [0.477209284166553];
% C_Hist_mean = [0.716885884283468];
% 
% Pd_Hist_std = [0.344955855390098];
% F_Hist_std = [0.00914224344405955];
% C_Hist_std = [0.00765435909684104];

% Plot K ratio w/ stds
fig = figure();
hold on;
x = 2000 : 10 : 2090;
%
neg = RCP45_Ratio - RCP45_Ratio_low;
pos = RCP45_Ratio_high - RCP45_Ratio;
errorbar(x, RCP45_Ratio, neg, pos, 'o-');
%
neg = RCP85_Ratio - RCP85_Ratio_low;
pos = RCP85_Ratio_high - RCP85_Ratio;
errorbar(x, RCP85_Ratio, neg, pos, 'o-');
%
xlabel('Decade');
ylabel('K Ratio');
legend('RCP 4.5', 'RCP 8.5');
%
exportgraphics(fig, ['Output_pd_f_c_calibration_decadal/K_Ratio_combined_errorbar.'  file_type], 'Resolution', file_res);
close fig 1;

% Plot horizontal spans
fig = figure();
hold on;
%
xi = 2010 : 10 : 2100;
x(1) = 2000;
j = 2;
for i = 1 : size(xi, 2)
    x(j) = xi(i);
    x(j + 1) = xi(i);
    j = j + 2;
end
x(end) = [];
%
posi = transpose(RCP45_Ratio_high); 
negi = transpose(RCP45_Ratio_low);
pos = NaN;
neg = NaN;
j = 1;
%
for i = 1 : size(xi, 2)
    pos(j) = posi(i);
    pos(j + 1) = posi(i);
    neg(j) = negi(i);
    neg(j + 1) = negi(i);
    j = j + 2;
end
%
x2 = [x, fliplr(x)];
inBetween = [pos, fliplr(neg)];
fill(x2, inBetween, 'b', 'FaceAlpha', 0.5);
%
posi = transpose(RCP85_Ratio_high); 
negi = transpose(RCP85_Ratio_low);
pos = NaN;
neg = NaN;
j = 1;
%
for i = 1 : size(xi, 2)
    pos(j) = posi(i);
    pos(j + 1) = posi(i);
    neg(j) = negi(i);
    neg(j + 1) = negi(i);
    j = j + 2;
end
%
x2 = [x, fliplr(x)];
inBetween = [pos, fliplr(neg)];
fill(x2, inBetween, 'r', 'FaceAlpha', 0.5);
%
meani = transpose(RCP45_Ratio); 
mean = NaN;
j = 1;
%
for i = 1 : size(xi, 2)
    mean(j) = meani(i);
    mean(j + 1) = meani(i);
    j = j + 2;
end
%
plot(x, mean, 'b');
%
meani = transpose(RCP85_Ratio); 
mean = NaN;
j = 1;
%
for i = 1 : size(xi, 2)
    mean(j) = meani(i);
    mean(j + 1) = meani(i);
    j = j + 2;
end
%
plot(x, mean, 'r');
%
plot([2000, 2100], [1, 1], 'k--');
%
xlabel('Year');
ylabel('K Ratio');
legend('RCP 4.5', 'RCP 8.5', 'location', 'northwest');
%
exportgraphics(fig, ['Output_pd_f_c_calibration_decadal/K_Ratio_combined_span.'  file_type], 'Resolution', file_res);
close fig 1;





%% More table exports

% Table
T45h = table(RCP45_Ratio_high);
T45l = table(RCP45_Ratio_low);
T85h = table(RCP85_Ratio_high);
T85l = table(RCP85_Ratio_low);
writetable(T45h, 'Output_pd_f_c_calibration_decadal/RCP45_High_Data.csv');
writetable(T45l, 'Output_pd_f_c_calibration_decadal/RCP45_Low_Data.csv');
writetable(T85h, 'Output_pd_f_c_calibration_decadal/RCP85_High_Data.csv');
writetable(T85l, 'Output_pd_f_c_calibration_decadal/RCP85_Low_Data.csv');
