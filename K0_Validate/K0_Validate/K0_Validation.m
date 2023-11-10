clear all
close all
clc




%% Parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dxy = [26.6736749983000]

m_sp = 0.503
n_sp = 1.224

min_area = 1  %10000000;  % minimum drainage area to prune network: m^2    WI=8

file_type = 'pdf'
file_res = 600

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get list
contents = dir( 'Created_Terrains/' );
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

% Convert list strings to numbers
list_convert = list;
list_num = [];
for i = 1 : size(list, 2);
    list_convert{i} = erase(list_convert{i}, '.tif');
    list_convert{i} = str2double(list_convert{i});
    list_num(1, i) = list_convert{i};
end

% Reorder
% list_copy{1, 1} = list{1, 2};
% list_copy{1, 2} = list{1, 3};
% list_copy{1, 3} = list{1, 1};
% list = list_copy;
list = flip(list)

% Make output directory
mkdir('Output');

% Set input directory
direct = 'Created_Terrains/';





%% Analysis real DEM

% Inport real DEM
%[DEMr,FDr,Ar,Sr] = MakeStreams('/Users/Chris/Dropbox/BC_Landlab/SPACE_2022/Chestatee/DEMs/Chestatee.tif', min_area, 'no_data_exp','auto');
[DEMr,FDr,Ar,Sr] = MakeStreams('C:/Users/sheehacz/Dropbox/BC_Landlab/SPACE_2022/New_Models/Input/Chestatee.tif', min_area, 'no_data_exp','auto');

% Slope-area of real DEM
fig = figure(1);
hold on;
SAr = slopearea(Sr, DEMr, Ar);
exportgraphics(fig, ['Output/SA_real.' file_type], 'Resolution', file_res);
close fig 1;

% Logrithm of SAr data
alogr = log(SAr.a);
glogr = log(SAr.g);

% Polyfit data
fit = fitlm(alogr, glogr);

% Transpose DEMr
zlistr = [];
for c = 1 : size(DEMr.Z, 2)
    zlistr = [zlistr; DEMr.Z(:, c)];
end

% Plot real DEM
fig = figure(1);
imageschs(DEMr, [], 'colormap', 'turbo');
exportgraphics(fig, ['Output/DEM_real.' file_type], 'Resolution', file_res);
close fig 1;

% Gamma distribution of elevations
pd = fitdist(zlistr, 'gamma');
y = pdf(pd, zlistr);

% Find elevation of max y
max(y);
find(y == ans);
zlistr(ans);
maxpdfr = mean(ans);

% Plot single pdf
fig = figure(1);
hold on;
plot(zlistr, y, '.');
xlabel('Elevation (m)');
ylabel('pdf');
exportgraphics(fig, ['Output/pdf_real.' file_type], 'Resolution', file_res);
close fig 1;

% Plot combined pdf
fig = figure(3);
hold on;
plot(zlistr, y, 'K-', 'LineWidth', 1);

% Calculate real Ksn
ksnr = ksn(Sr, DEMr, Ar, m_sp / n_sp);
ksnr = smooth(Sr, ksnr);

% Map real Ksn values
fig = figure(1);
hold on;
imageschs(DEMr, [], 'colorbar', false, 'colormap', 'gray');
scatter(Sr.x, Sr.y, 5, ksnr, 'filled'); 
caxis([0 max(ksnr)]);
colorbar;
exportgraphics(fig, ['Output/ksn_real.' file_type], 'Resolution', file_res);
close fig 1;

% Record values
mr = fit.Coefficients{2, 1};
br = fit.Coefficients{1, 1};
r2r = fit.Rsquared.Adjusted;
thetar = SAr.theta;
ksr = SAr.ks; 
zmeanr = nanmean(nanmean(DEMr.Z));
zstdr = nanstd(nanstd(DEMr.Z));
zmaxr = nanmax(nanmax(DEMr.Z));
zmedianr = nanmedian(nanmedian(DEMr.Z));

% Create empty arrays
m = ones(size(list, 2)) * NaN;
b = ones(size(list, 2)) * NaN;
r2 = ones(size(list, 2)) * NaN;
theta = ones(size(list, 2)) * NaN;
ks = ones(size(list, 2)) * NaN;
zmean = ones(size(list, 2)) * NaN;
zstd = ones(size(list, 2)) * NaN;
zmax = ones(size(list, 2)) * NaN;
zmedian = ones(size(list, 2)) * NaN;

% % Calculate real chi
% fig = chiplot(Sr, DEMr, Ar, 'a0', min_area, 'mn', m_sp / n_sp);
% exportgraphics(figure(1), ['Output/chi_real.' file_type], 'Resolution', file_res);
% close fig 1;

% Plot combined SA
fig = figure(4);
hold on
p = polyfit(log10(SAr.a), log10(SAr.g), 1);
y = polyval(p, log10(SAr.a));
hold on;
plot(10.^(log10(SAr.a)), 10.^(y), 'DisplayName', ['Real'], 'LineWidth', 1, 'Color', 'k');
plot(10.^(log10(SAr.a)), 10.^(log10(SAr.g)), 's', 'Color', 'k');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log')

% Boxplot
B(1:size(zlistr, 1), 1) = zlistr;





%% Analysis model DEMs

for i = 1 : size(list, 2)
    
    % Makestreams
    DEMm = GRIDobj([direct list{i}]);
    
    % Handle no-data values
    for r = 1 : size(DEMm.Z, 1)
        for c = 1 : size(DEMm.Z, 2)
            if DEMm.Z(r,c) == -99999
                DEMm.Z(r,c) = NaN;
            end
        end
    end
    
    % Create FD, A, and S for each DEM
    DEMm = fillsinks(DEMm)
    FDm = FLOWobj(DEMm);
    Am = flowacc(FDm);
    Sm = STREAMobj(FDm, Am >= min_area); 
    
    % Slope-area relationship of model
    fig = figure(1);
    SAm = slopearea(Sm, DEMm, Am);
    exportgraphics(fig, ['Output/SA_' list{i} '.' file_type], 'Resolution', file_res)
    close fig 1;
    
    % Logrithm of SAr data
    alogm = log(SAm.a);
    glogm = log(SAm.g);

    % Polyfit data
    fit = fitlm(alogm, glogm);
    
    % Record values
    m(1, i) = fit.Coefficients{2, 1};
    b(1, i) = fit.Coefficients{1, 1};
    r2(1, i) = fit.Rsquared.Adjusted;
    theta(1, i) = SAm.theta;
    ks(1, i) = SAm.ks;
    zmean(1, i) = nanmean(nanmean(DEMm.Z));
    zstd(1, i) = nanstd(nanstd(DEMm.Z));
    zmax(1, i) = nanmax(nanmax(DEMm.Z));
    zmedian(1, i) = nanmedian(nanmedian(DEMm.Z));
    
    % Transpose DEMm
    zlistm = [];
    for c = 1 : size(DEMm.Z, 2)
        zlistm = [zlistm; DEMm.Z(:, c)];
    end  
    
    % Plot model DEM
    fig = figure(1);
    imageschs(DEMm, [], 'colormap', 'turbo');
    exportgraphics(fig, ['Output/DEM_' list{i} '.' file_type], 'Resolution', file_res);
    close fig 1;  
    
    % Plot model vs real DEM
    fig = figure(1);
    subplot(1, 2, 1);
    hold on;
    imageschs(DEMr, [], 'colormap', 'turbo');
    %caxis([min(zlistr) max(zlistr)]);
    colorbar;
    xlabel('Real')
    subplot(1, 2, 2);
    hold on;
    imageschs(DEMm, [], 'colormap', 'turbo');
    %caxis([min(zlistr) max(zlistr)]);
    colorbar;
    xlabel('Model')
    exportgraphics(fig, ['Output/DEM_comparison_' list{i} '.' file_type], 'Resolution', file_res);
    close fig 1;         
    
    % Gamma distribution
    pd = fitdist(zlistm, 'gamma');
    y = pdf(pd, zlistm);
    
    % Find elevation of max y
    max(y);
    find(y == ans);
    zlistm(ans);
    maxpdf(1, i) = mean(ans);
    
    % Plot single
    fig = figure(1);
    hold on;
    plot(zlistm, y, '.');
    xlabel('Elevation (m)');
    ylabel('pdf');
    %title(['K = ' list{i}]);
    exportgraphics(fig, ['Output/pdf_' list{i} '.' file_type], 'Resolution', file_res);
    close fig 1;

    % Plot Combined
    fig = figure(3);
    hold on;
    % plot(zlistm, y, '.', 'MarkerSize', 1, 'DisplayName', ['K = ' erase(list{i}, '.tif')]);
    plot(zlistm, y, '-', 'LineWidth', 1);

    % Calculate model Ksn
    ksnm = ksn(Sm ,DEMm ,Am , m_sp / n_sp);
    ksnm = smooth(Sm, ksnm);

%     % Map model Ksn values
%     fig = figure(1);
%     hold on;
%     imageschs(DEMm, [], 'colorbar', false, 'colormap', 'gray');
%     scatter(Sm.x, Sm.y, 5, ksnm, 'filled'); 
%     caxis([0 max(ksnr)]);
%     colorbar;
%     exportgraphics(fig, ['Output/ksn_' list{i} '.' file_type], 'Resolution', file_res);
%     close fig 1;
    
    % Plot combined SA
    fig = figure(4);
    hold on
    p = polyfit(log10(SAm.a), log10(SAm.g), 1);
    y = polyval(p, log10(SAm.a));
    hold on;
    plot(SAm.a, 10.^(y), 'DisplayName', ['K = ' list{i}]);
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log')
    %plot(log(SAm.a), log(SAm.g), 's');
    
    % Boxplot
    B(1:size(zlistm, 1), i + 1) = zlistm;

    
end

% Finish combined pdf
%uistack(h,'top');
fig = figure(3);
hold on;
xlabel('Elevation (m)');
ylabel('pdf');
%legend();
legend('Real landscape', 'Model Mean K0', 'Model Min K0', 'Model Max K0');
exportgraphics(fig, ['Output/pdf_combined.' file_type], 'Resolution', file_res);
close fig 3;

% Finish combined SA
fig = figure(4);
xlabel('log( drainage area )')
ylabel('log( slope )')
legend();
exportgraphics(fig, ['Output/SA_combined.' file_type], 'Resolution', file_res)
close fig 4;

% Finish Boxplot
fig = figure(5)
boxplot(B, 'Labels', {'DEM', 'Model K0', 'Model Min K0', 'Model Max K0'});
set(gca,'XTickLabel',{'DEM', 'Model K0', 'Model Min K0', 'Model Max K0'});
exportgraphics(fig, ['Output/Boxplot.' file_type], 'Resolution', file_res)
close fig 5;

% Finish Boxchart
fig = figure(6)
boxchart(B, 'MarkerStyle','none');
set(gca,'XTickLabel',{'Real landscape', 'K_0', 'K_0 min', 'K_0 max'}, 'fontsize', 11);
ylim([300, 800])
ylabel('Elevation (m)')
exportgraphics(fig, ['Output/Boxchart.' file_type], 'Resolution', file_res)
close fig 6;



%% Final Plots

% Remove extra values
m(2 : end, :) = [];
b(2 : end, :) = [];
r2(2 : end, :) = [];
theta(2 : end, :) = [];
ks(2 : end, :) = [];
zmean(2 : end, :) = [];
zstd(2 : end, :) = [];
zmax(2 : end, :) = [];
zmedian(2 : end, :) = [];

% Sort list_num
[sorted, indices] = sort(list_num);

% % Plot m
% fig = figure(1);
% hold on;
% semilogx(list_num(indices), m(indices), '.-');
% y = ones(length(list_num)) * mr;
% y(2, :) = [];
% semilogx(list_num, y);
% xlabel('K');
% ylabel('m')
% legend('model', 'real');
% exportgraphics(fig, ['Output/m.' file_type], 'Resolution', file_res);
% close fig 1;
% 
% % Plot b
% fig = figure(1);
% hold on;
% plot(list_num(indices), b(indices), '.-');
% y = ones(length(list_num)) * br;
% y(2, :) = [];
% plot(list_num, y);
% xlabel('K');
% ylabel('b')
% legend('model', 'real');
% 
% exportgraphics(fig, ['Output/b.' file_type], 'Resolution', file_res);
% close fig 1;
% 
% % Plot r2
% fig = figure(1);
% hold on;
% plot(list_num(indices), r2(indices), '.-');
% y = ones(length(list_num)) * r2r;
% y(2, :) = [];
% plot(list_num, y);
% xlabel('K');
% ylabel('r-squared')
% legend('model', 'real');
% exportgraphics(fig, ['Output/r2.' file_type], 'Resolution', file_res);
% close fig 1;
% 
% % Plot theta
% fig = figure(1);
% hold on;
% plot(list_num(indices), theta(indices), '.-');
% y = ones(length(list_num)) * thetar;
% y(2, :) = [];
% plot(list_num, y);
% xlabel('K');
% ylabel('theta')
% legend('model', 'real');
% exportgraphics(fig, ['Output/theta.' file_type], 'Resolution', file_res);
% close fig 1;
% 
% % Plot ks
% fig = figure(1);
% hold on;
% plot(list_num(indices), ks(indices), '.-');
% y = ones(length(list_num)) * ksr;
% y(2, :) = [];
% plot(list_num, y);
% xlabel('K');
% ylabel('ks')
% legend('model', 'real');
% exportgraphics(fig, ['Output/ks.' file_type], 'Resolution', file_res);
% close fig 1;

% Fit
p = polyfit(list_num(indices), zmean(indices), 2);
a = p(1);
b = p(2);
c = p(3) - zmeanr;
sqrt = ((b ^ 2) - (4 * a * c)) ^ (0.5);
K_elev_fit = (-b - sqrt) / (2 * a);

% % Plot mean elev
% fig = figure(1);
% hold on;
% x2 = [list_num(indices), fliplr(list_num(indices))];
% curve1 = zmean * NaN;
% curve1(1, :) = zmeanr - zstdr;
% curve2 = curve1 * NaN;
% curve2(1, :) = zmeanr + zstdr;
% shade = [curve1, fliplr(curve2)];
% fill(x2, shade, 'k', 'FaceAlpha', 0.25);
% errorbar(list_num(indices), zmean(indices), zstd(indices), '*-'); % 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k'
% plot(x2, ones(1, size(x2, 2)) * zmeanr, 'k');
% xline(K_elev_fit, 'k--');
% xlabel('K');
% ylabel('Mean elevation')
% legend('Real mean elevation (+- 1 std)', 'Modeled mean elevation (+- 1 std)', '', 'Best-fit K');
% exportgraphics(fig, ['Output/mean_elevation.' file_type], 'Resolution', file_res);
% close fig 1;

% % Plot max elev
% fig = figure(1);
% hold on;
% plot(list_num(indices), zmax(indices), '-s', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k');
% y = ones(length(list_num)) * zmaxr;
% y(2, :) = [];
% plot(list_num, y);
% xlabel('K');
% ylabel('Max elevation')
% legend('model', 'real');
% exportgraphics(fig, ['Output/max_elevation.' file_type], 'Resolution', file_res);
% close fig 1;

% % Plot median elev
% fig = figure(1);
% hold on;
% plot(list_num(indices), zmedian(indices), '-s', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k');
% y = ones(length(list_num)) * zmedianr;
% y(2, :) = [];
% plot(list_num, y);
% xlabel('K');
% ylabel('Median elevation')
% legend('model', 'real');
% exportgraphics(fig, ['Output/median_elevation.' file_type], 'Resolution', file_res);
% close fig 1;

% Calculate best fit K (based on interpolation of curve in "mean_elevation.file_type")




%%

% % Mean minimization
% minimization_mean = abs(zmeanr - zmean);
% fig = figure(1);
% hold on;
% plot(list_num(indices), minimization_mean(indices), '-s', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k');
% xlabel('K');
% ylabel('Mean Elevation Minimization')
% exportgraphics(fig, ['Output/minimization_mean.' file_type], 'Resolution', file_res);
% close fig 1;
% 
% % Median minimization
% minimization_median = abs(zmedianr - zmedian);
% fig = figure(1);
% hold on;
% plot(list_num(indices), minimization_median(indices), '-s', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k');
% xlabel('K');
% ylabel('Median Elevation Minimization')
% exportgraphics(fig, ['Output/minimization_median.' file_type], 'Resolution', file_res);
% close fig 1;
% 
% % Pdf minimization
% minimization_pdf = abs(maxpdfr - maxpdf);
% fig = figure(1);
% hold on;
% plot(list_num(indices), minimization_pdf(indices), '-s', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k');
% xlabel('K');
% ylabel('PDF Minimization')
% exportgraphics(fig, ['Output/minimization_pdf.' file_type], 'Resolution', file_res);
% close fig 1;





%% For CSDMS Poster

% % Makestreams
% DEMm = GRIDobj([direct '7.5E-7.tif']);
% 
% % Handle no-data values
% for r = 1 : size(DEMm.Z, 1)
%     for c = 1 : size(DEMm.Z, 2)
%         if DEMm.Z(r,c) == -99999
%             DEMm.Z(r,c) = NaN;
%         end
%     end
% end
% 
% % Create FD, A, and S for each DEM
% DEMm = fillsinks(DEMm)
% FDm = FLOWobj(DEMm);
% Am = flowacc(FDm);
% Sm = STREAMobj(FDm, Am >= min_area);
% 
% % Plot model vs real DEM
% fig = figure(1);
% subplot(2, 1, 1);
% hold on;
% imageschs(DEMr, [], 'colormap', 'turbo', 'nancolor', [0 0 0]);
% %caxis([min(zlistr) max(zlistr)]);
% colorbar;
% xlabel('Real')
% subplot(2, 1, 2);
% hold on;
% imageschs(DEMm, [], 'colormap', 'turbo', 'nancolor', [0 0 0]);
% %caxis([min(zlistr) max(zlistr)]);
% colorbar;
% xlabel('Model')
% exportgraphics(fig, ['Output/CSDMS_A' list{i} '.' file_type], 'Resolution', file_res);
% close fig 1;         















