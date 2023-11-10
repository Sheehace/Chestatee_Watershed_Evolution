clear all
close all
clc




%% Parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dxy = [26.6736749983000]

%Basin = 'A'

%m_sp = 0.503 
%n_sp = 1.224

min_area = 1   %10000000;  % minimum drainage area to prune network: m^2    WI=8

file_type = 'png'
file_res = 900

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get list
contents = dir( ['Created_Terrains/'] );
names = {contents.name};
list = {};
i = 1;
for c = 1 : size(names, 2)
    if contains(names{c}, '.tif') == 1 
        if contains(names{c}, '.aux') == 0 
            list{i} = names{c};
            i = i + 1;
        end
    end
end

% Delete "A.tif" from list
%list(end) = [];

% Convert list strings to numbers
list_convert = list;
list_num = [];
for i = 1 : size(list, 2);
    list_convert{i} = erase(list_convert{i}, '.tif');
    list_convert{i} = str2double(list_convert{i});
    list_num(1, i) = list_convert{i};
end

% Make output directory
mkdir('Output');

% Set input directory
direct = 'Created_Terrains/';





%% Analysis

% Inport Chestatee DEM
%[DEMc,FDc,Ac,Sc] = MakeStreams('/Users/Chris/Dropbox/BC_Landlab/SPACE_2022/New_Models/Input/Chestatee.tif', min_area, 'no_data_exp','auto');
%[DEMc,FDc,Ac,Sc] = MakeStreams('/Users/Chris/Dropbox/BC_Landlab/AP_Scarp/Input/YR.tif', min_area, 'no_data_exp','auto');
[DEMc,FDc,Ac,Sc] = MakeStreams('C:/Users/sheehacz/Dropbox/BC_Landlab/SPACE_2022/New_Models/Input/Chestatee.tif', min_area, 'no_data_exp','auto');

% Based on inspection of Chestatee DEM
min_slope = 0.04;
min_grad = (dxy^2) * 2;

% Slope-area of Chestatee DEM
fig = figure(1);
hold on;
SAc = slopearea(Sc, DEMc, Ac);
exportgraphics(fig, ['Output/SA_Chestatee.' file_type], 'Resolution', file_res);
close fig 1;

% Plot Chestatee DEM
fig = figure(1);
imageschs(DEMc, [], 'colormap', 'turbo');
exportgraphics(fig, ['Output/DEM_Chestatee.' file_type], 'Resolution', file_res);
close fig 1;

% Find values to use in regression
use1 = SAc.g > min_slope;
use2 = SAc.a > min_grad;
use = use1 & use2;

% Regress slope-area
xc = log(SAc.a(use));
yc = log(SAc.g(use));
pc = polyfit(xc, yc, 1);





%% Loop through DEMs

% Sort list_num
[sorted, indices] = sort(list_num);

% Initiate variables
Rm_mean = [];
Rm_std = [];
maxg_index = [];

for i = 1 : size(list, 2)
    
    % Makestreams
    DEMi = GRIDobj([direct list{i}]);
    
    % Handle no-data values
    for r = 1 : size(DEMi.Z, 1)
        for c = 1 : size(DEMi.Z, 2)
            if DEMi.Z(r,c) == -99999
                DEMi.Z(r,c) = NaN;
            end
        end
    end
    
    % Create FD, A, and S for each DEM
    DEMi = fillsinks(DEMi)
    FDm = FLOWobj(DEMi);
    Am = flowacc(FDm);
    Sm = STREAMobj(FDm, Am >= (min_area / (dxy^2))); 
    
    % Slope-area relationship of model
    fig = figure(1);
    SAi = slopearea(Sm, DEMi, Am);
    exportgraphics(fig, ['Output/SA_' list{i} '.' file_type], 'Resolution', file_res)
    close fig 1;
    
    % Regress slope-area
    xi = log(SAi.a(use));
    yi = log(SAi.g(use));
    pi = polyfit(xi, yi, 1);
    
    % Plot model DEM
    fig = figure(1);
    imageschs(DEMi, [], 'colormap', 'turbo');
    exportgraphics(fig, ['Output/DEM_' list{i} '.' file_type], 'Resolution', file_res);
    close fig 1;  
    
    % Plot combined DEMs
    fig = figure(1);
    hold on;
    subplot(1, 2, 1)
    imageschs(DEMc, [], 'colormap', 'turbo');
    xlabel('real');
    subplot(1, 2, 2)
    hold on;
    imageschs(DEMi, [], 'colormap', 'turbo');
    xlabel('model');
    exportgraphics(fig, ['Output/DEM_Comparison_' list{i} '.' file_type], 'Resolution', file_res);
    close fig 1;  
    
    % Plot combined SAs
    fig = figure();
    hold on;
    plot(log(SAc.a), log(SAc.g), 'ks');
    plot(log(SAc.a), (log(SAc.a) * pc(1)) + pc(2), 'k');
    plot(log(SAi.a), log(SAi.g), 'bs');
    plot(log(SAi.a), (log(SAi.a) * pi(1)) + pi(2), 'b');
    exportgraphics(fig, ['Output/SA_Comparison_' list{i} '.' file_type], 'Resolution', file_res);
    close fig 1;  

    % Plot combined SAs, no trendline
    fig = figure();
    hold on;
    plot(log(SAc.a), log(SAc.g), 'ks');
    % plot(log(SAc.a), (log(SAc.a) * pc(1)) + pc(2), 'k');
    plot(log(SAi.a), log(SAi.g), 'bs');
    % plot(log(SAi.a), (log(SAi.a) * pi(1)) + pi(2), 'b');
    exportgraphics(fig, ['Output/SA_Comparison_clean_' list{i} '.' file_type], 'Resolution', file_res);
    close fig 1;  
    
    % Find index of max g in SA object
    maxg_index(1, i) = find(SAi.g == max(SAi.g));
    if SAi.g(1) > SAi.g(3);
        j = 1;
    elseif SAi.g(1) < SAi.g(3);
        j = 3;
    end
    one_or_three_larger_slope(1, i) = j;

    % Store data to create figure for paper
    if string(list{i}) == '200.tif'
        X_low = SAi.a;
        Y_low = SAi.g;
    end
    %
    if string(list{i}) == '325.tif'
        X_high = SAi.a;
        Y_high = SAi.g;
    end


end

In_Range = maxg_index == 2 & one_or_three_larger_slope == 1;

% Create Figure for paper
fig = figure()
hold on;
plot(SAc.a, SAc.g, 'ks');
plot(X_low, Y_low, 'bs');
plot(X_high, Y_high, 'rs');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('Drainage area (m^2)');
ylabel('Slope');
legend('Real', 'D* = 200', 'D* = 325');
exportgraphics(fig, ['Output/For_Paper.' file_type], 'Resolution', file_res);
close fig 1;  


%%

% % Plot roughness and diff
% Rm_diff = abs(Rm_mean - Ra_mean);
% %
% fig = figure(1);
% subplot(2, 1, 1);
% hold on;
% %
% x2 = [list_num(indices), fliplr(list_num(indices))];
% curve1 = Rm_mean * NaN;
% curve1(1, :) = 0 - Ra_std;
% curve2 = Rm_mean * NaN;
% curve2(1, :) = 0 + Ra_std;
% shade = [curve1, fliplr(curve2)];
% fill(x2, shade, 'g');
% %
% errorbar(list_num(indices), Rm_diff(indices), Rm_std(indices), '-s', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')
% xlabel('D/K');
% ylabel('abs[ (real ri) - (model ri)');
% 
% %
% fig = figure(1);
% subplot(2, 1, 2);
% hold on;
% %
% x2 = [list_num(indices), fliplr(list_num(indices))];
% curve1 = Rm_mean * NaN;
% curve1(1, :) = Ra_mean - Ra_std;
% curve2 = Rm_mean * NaN;
% curve2(1, :) = Ra_mean + Ra_std;
% shade = [curve1, fliplr(curve2)];
% fill(x2, shade, 'g');
% %
% errorbar(list_num(indices), Rm_mean(indices), Rm_std(indices), '-s', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k')
% xlabel('D/K');
% ylabel('roughness index');
% %
% exportgraphics(fig, ['Output/Roughness.' file_type], 'Resolution', file_res);
% close fig 1;





%% For CSDMS Poster
% 
% fig = figure(1);
% %
% subplot(2, 3, 5);
% hold on;
% DEMm = GRIDobj([direct '80.tif']);
% DEMm = fillsinks(DEMm)
% for r = 1 : size(DEMm.Z, 1)
%     for c = 1 : size(DEMm.Z, 2)
%         if DEMm.Z(r,c) == -99999
%             DEMm.Z(r,c) = NaN;
%         end
%     end
% end
% imageschs(DEMm, [], 'colormap', 'turbo', 'nancolor', [0 0 0]);
% title('D/K = 80');
% %
% subplot(2, 3, 4);
% hold on;
% DEMm = GRIDobj([direct '800.tif']);
% DEMm = fillsinks(DEMm)
% for r = 1 : size(DEMm.Z, 1)
%     for c = 1 : size(DEMm.Z, 2)
%         if DEMm.Z(r,c) == -99999
%             DEMm.Z(r,c) = NaN;
%         end
%     end
% end
% imageschs(DEMm, [], 'colormap', 'turbo', 'nancolor', [0 0 0]);
% title('D/K = 800 (optimized)');
% %
% subplot(2, 3, 6);
% hold on;
% DEMm = GRIDobj([direct '8000.tif']);
% DEMm = fillsinks(DEMm)
% for r = 1 : size(DEMm.Z, 1)
%     for c = 1 : size(DEMm.Z, 2)
%         if DEMm.Z(r,c) == -99999
%             DEMm.Z(r,c) = NaN;
%         end
%     end
% end
% imageschs(DEMm, [], 'colormap', 'turbo', 'nancolor', [0 0 0]);
% title('D/K = 8000');
% %
% subplot(2, 3, 1);
% imageschs(DEMa, [], 'colormap', 'turbo', 'nancolor', [0 0 0]);
% title('Real')
% %
% exportgraphics(fig, ['Output/CSDMS_A.' file_type], 'Resolution', file_res);
% close fig 1;
% 
% 
% 
% 
% 
% 
% 
% 
% DEMm = GRIDobj([direct '800.tif']);
% for r = 1 : size(DEMm.Z, 1)
%     for c = 1 : size(DEMm.Z, 2)
%         if DEMm.Z(r,c) == -99999
%             DEMm.Z(r,c) = NaN;
%         end
%     end
% end
% DEMm = fillsinks(DEMm)
% FDm = FLOWobj(DEMm);
% Am = flowacc(FDm);
% Sm = STREAMobj(FDm, Am >= (min_area / (dxy^2))); 
% %
% fig = figure(1);
% %
% subplot(1, 2, 1);
% hold on;
% slopearea(Sc, DEMc, Ac);
% title('Real');
% xlabel('drainage area (m2)');
% %
% subplot(1, 2, 2);
% hold on;
% slopearea(Sm, DEMm, Am);
% title('D/K = 800 (optimized)');
% xlabel('drainage area (m2)');
% %
% set(gcf, 'units','inches', 'position', [0,0, 8, 4])
% exportgraphics(fig, ['Output/CSDMS_B.' file_type], 'Resolution', file_res);
% close fig 1;
























