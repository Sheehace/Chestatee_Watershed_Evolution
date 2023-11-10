#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 15:19:30 2022

@author: Chris
"""
# %% Initiate

# ENTER VARIABLES ############################################################
##############################################################################

import inspect
import os
Directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
project_directory = Directory.replace('\k1_validate', '')

from landlab.io.esri_ascii import read_esri_ascii
from landlab.plot.graph import plot_graph
from matplotlib import cm
from matplotlib import pyplot as plt
from landlab.io.esri_ascii import write_esri_ascii
from landlab.components import StreamPowerEroder, Space, SpaceLargeScaleEroder, DepressionFinderAndRouter, LinearDiffuser, TaylorNonLinearDiffuser,  FlowAccumulator, ChannelProfiler, SteepnessFinder, ChiFinder, Lithology, LithoLayers, NormalFault
from landlab.plot.video_out import VideoPlotter
from landlab import RasterModelGrid, imshow_grid
import sys
from os import path
import os.path
import os
import math
import pandas as pd
from numpy import nan, isnan
import numpy as np
#Directory = 'C:/Users/sheehacz/Manuscript_Codes/K1_Validate'

# Gauge_Hydro_path = 'C:/Users/sheehacz/Manuscript_Codes/Klc_Validate/02333500_discharge.csv'
basin_area = 608568641.6459954      # m^2, entire Chestaee
# m^2, calculated from 153 mi^2 reported on USGS Waterwatch # 246229632
gauge_upstream_area = 396268179.3

# Cum_02333500.csv
# Gauge_Sed_path = 'C:/Users/sheehacz/Manuscript_Codes/K1_Validate/02333500_USE.csv'
Gauge_Sed_path = Directory + '/02333500_USE.csv'



Export_format = 'pdf'
dpi = 900

e = 2.71828

# ENTER VARIABLES ############################################################
##############################################################################

# Print space
print(' ')

# Import libraries
print('Importing libraries...')
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
warnings.filterwarnings("ignore")

# Create directory
print('Creating directory...')
if path.exists(str(Directory)+'/Output') == False:
    os.mkdir(str(Directory)+'/Output')

# Import and clean gauge hydro data
print('Importing and cleaning Gauge_Hydro...')
skiprows = np.arange(0, 36)
skiprows = np.append(skiprows, 37)
# Gauge_Hydro = pd.read_csv(Gauge_Hydro_path, skiprows=skiprows)

# # Convert gauge discharge (cfs to m^3 per year) and height (ft to m)
# Gauge_Hydro["daily_discharge"] *= 893000.074
# Gauge_Hydro["daily_mean_gauge_height"] *= 0.3048

# If using Cum_02333500.csv
# # Import and clean Gauge_Sed data
# print('Importing and cleaning Gauge_Sed data...')
# Gauge_Sed = pd.read_csv(Gauge_Sed_path, skiprows = [1, 2]) # nrows = 545
# for i in Gauge_Sed.columns:
#     for j in np.arange(0, len(Gauge_Sed)):
#         if type(Gauge_Sed[i][j]) == str:
#             if "E" in Gauge_Sed[i][j]:
#                 Gauge_Sed[i][j] = nan
#             if Gauge_Sed[i][j] == "#DIV/0!":
#                 Gauge_Sed[i][j] = nan
#             if Gauge_Sed[i][j] == "< 1":
#                 Gauge_Sed[i][j] = nan
#             if Gauge_Sed[i][j] == "< 1 ":
#                 Gauge_Sed[i][j] = nan
#         if (type(Gauge_Sed[i][j]) == str):
#             Gauge_Sed[i][j] = float(Gauge_Sed[i][j])
# Gauge_Sed['sediment_flux_m3_year'] = Gauge_Sed['Suspended sediment discharge, short tons per day'] * (907.18474) * (1 / 2650) * (365)   # Convert short tons per day to cubic meters per year. 907.18474 kg per short ton, 2650 kg per cubic meter (quartz), 365 days per year
# Gauge_Sed['water_discharge_m3_year'] = Gauge_Sed['Enforced Discharge, cfs'] * 893000.074    # Convert cfs to m^3 / yr.

# If using 02333500_USE.csv (Reccommended)
# Import and clean Gauge_Sed data
print('Importing and cleaning Gauge_Sed data...')
Gauge_Sed = pd.read_csv(Gauge_Sed_path, skiprows=[1, 2])  # nrows = 545
# for i in Gauge_Sed.columns:
#     for j in np.arange(0, len(Gauge_Sed)):
#         if type(Gauge_Sed[i][j]) == str:
#             if "E" in Gauge_Sed[i][j]:
#                 Gauge_Sed[i][j] = nan
#             if Gauge_Sed[i][j] == "#DIV/0!":
#                 Gauge_Sed[i][j] = nan
#             if Gauge_Sed[i][j] == "< 1":
#                 Gauge_Sed[i][j] = nan
#             if Gauge_Sed[i][j] == "< 1 ":
#                 Gauge_Sed[i][j] = nan
#         if (type(Gauge_Sed[i][j]) == str):
#             Gauge_Sed[i][j] = float(Gauge_Sed[i][j])
# Convert short tons per day to cubic meters per year. 907.18474 kg per short ton, 2650 kg per cubic meter (quartz), 365 days per year
Gauge_Sed['sediment_flux_m3_year'] = Gauge_Sed['Suspended sediment discharge, short tons per day'] * \
    (907.18474) * (1 / 2650) * (365)
# Convert cfs to m^3 / yr.
Gauge_Sed['water_discharge_m3_year'] = Gauge_Sed['Enforced Discharge, cfs'] * 893000.074


# 10Be data
# SAP56_mean = 9.5 * 1E-6
# SAP56_std = 0.75 * 1E-6
# SAP56_recalc = 7.96 * 1E-6
# reusser_mean = np.array([11.38, 12.44, 17.16]) * 1E-6
# reusser_std = np.array([0.96, 0.95, 1.32]) * 1E-6
# reusser_recalc = np.array([9.15, 10.05]) * 1E-6

SAP56_mean_recalc = 7.96 * 1E-6          # m / yr^-1
SAP56_std_recalc = 0.63 * 1E-6        # m / yr^-1
SAP56_mean = 9.5 * 1E-6
SAP56_std = 0.75 * 1E-6
SAP60_mean = 17.16 * 1E-6
SAP60_std = 1.32 * 1E-6
Trimble_yield = 5.7 * 1E-5        # mm/100 years to m / yr
Trimble_upland = 120 * 1E-5         # # mm/100 years to m / yr

# Hydro = pd.read_csv('C:/Users/sheehacz/Manuscript_Codes/K1_Validate/02333500_Large.csv')
Hydro = pd.read_csv(Directory + '/02333500_Large.csv')

Hydro['daily_discharge_m3yr-1'] = nan
Hydro['daily_discharge_m3yr-1'] = Hydro['36255_00060_00003'] * \
    893000.074          # cfs to m3yr-1
# Hydro['year'] = nan
# year = 30
# for i in np.arange(1, np.size(Hydro['datetime'])):
#     substring = '/' + str(year)
#     Hydro['year'][i] = year + 1900


#     if Hydro['datetime'][i] == '12/31' + substring:
#         year += 1
#     print(i / np.size(Hydro['datetime']))


# %% Create rating curve

# Extract log values of gauge sed data
Qw = Gauge_Sed['water_discharge_m3_year'].values
Qw = Qw.astype(float)
Qw = np.log(Qw)
#
Qs = Gauge_Sed['sediment_flux_m3_year'].values
Qs = Qs.astype(float)
Qs = np.log(Qs)

# Remove nans
delete = np.ones(np.size(Qw)) * nan
for i in np.arange(0, np.size(Qw)):
    if isnan(Qw[i]) == True or isnan(Qs[i]) == True:
        delete[i] = True
    else:
        delete[i] = False
delete = np.where(delete == True)
delete = delete[0]

Qw = np.delete(Qw, delete)
Qs = np.delete(Qs, delete)

# Fit Regression.                           # log(Qw) = ( m * log(Qs) ) + b
polyfit = np.polyfit(Qw, Qs, deg=1)
m = polyfit[0]
b = polyfit[1]

# Create xaxis values for regression plot
line_x = np.array([np.min(Qw), np.max(Qw)])
line_y = (line_x * m) + b

# Plot
if path.exists(str(Directory)+'/Output') == False:
    os.mkdir(str(Directory)+'/Output')
plt.ioff()
fig = plt.figure()
#plt.plot(Qw, Qs, '.', zorder = 1)
plt.scatter(Qw, Qs, 5, zorder=1, c=Gauge_Sed['Year'])
plt.plot(line_x, line_y)
# plt.xlabel('$log(Qw)$ $(m^{{3}}yr^{-1})$')
# plt.ylabel('$log(Qs)$ $(m^{3}yr^{-1})$')
plt.xlabel('log(water discharge) $(m^{{3}}$ $yr^{-1})$')
plt.ylabel('log(sediment discharge) $(m^{3}$ $yr^{-1})$')
plt.text(20.5, 4, 'y = ' + str(m) + 'x ' + str(b))
plt.grid()
plt.tight_layout()
fig.savefig(str(Directory)+'/Output/Rating_Curve.' + Export_format, dpi=dpi)
plt.close(fig)

# Use Rating Curve
if path.exists(str(Directory)+'/Output') == False:
    os.mkdir(str(Directory)+'/Output')
plt.ioff()
fig = plt.figure()
#plt.plot(Qw, Qs, '.', zorder = 1)
plt.scatter(Qw, Qs, 5, zorder=1)
plt.plot(line_x, line_y, 'k--')
# plt.xlabel('$log(Qw)$ $(m^{{3}}yr^{-1})$')
# plt.ylabel('$log(Qs)$ $(m^{3}yr^{-1})$')
plt.xlabel('log(water discharge) $(m^{{3}}$ $yr^{-1})$')
plt.ylabel('log(sediment discharge) $(m^{3}$ $yr^{-1})$')
# plt.text(20.5, 4, 'y = ' + str(m) + 'x ' + str(b))
plt.grid()
plt.tight_layout()
fig.savefig(str(Directory)+'/Output/Rating_Curve_Use.' + Export_format, dpi=dpi)
plt.close(fig)

# Plot in groups
if path.exists(str(Directory)+'/Output') == False:
    os.mkdir(str(Directory)+'/Output')
plt.ioff()
fig = plt.figure()

plt.scatter(Qw[0: 21], Qs[0: 21], 5, zorder=1, color='blue')
polyfit = np.polyfit(Qw[0: 21], Qs[0: 21], deg=1)
m = polyfit[0]
b = polyfit[1]
line_x = np.array([np.min(Qw[0: 21]), np.max(Qw[0: 21])])
line_y = (line_x * m) + b
plt.plot(line_x, line_y, color='blue')

plt.scatter(Qw[21: 59], Qs[21: 59], 5, zorder=1, color='green')
polyfit = np.polyfit(Qw[21: 59], Qs[21: 59], deg=1)
m = polyfit[0]
b = polyfit[1]
line_x = np.array([np.min(Qw[21: 59]), np.max(Qw[21: 59])])
line_y = (line_x * m) + b
plt.plot(line_x, line_y, color='green')

plt.scatter(Qw[59: 159], Qs[59: 159], 5, zorder=1, color='orange')
polyfit = np.polyfit(Qw[59: 159], Qs[59: 159], deg=1)
m = polyfit[0]
b = polyfit[1]
line_x = np.array([np.min(Qw[59: 159]), np.max(Qw[59: 159])])
line_y = (line_x * m) + b
plt.plot(line_x, line_y, color='orange')

plt.scatter(Qw[159: 273], Qs[159: 273], 5, zorder=1, color='red')
polyfit = np.polyfit(Qw[159: 273], Qs[159: 273], deg=1)
m = polyfit[0]
b = polyfit[1]
line_x = np.array([np.min(Qw[159: 273]), np.max(Qw[159: 273])])
line_y = (line_x * m) + b
plt.plot(line_x, line_y, color='red')
plt.legend(['1957 - 1963 (n = 21)', '', '1972 - 1976 (n = 38)', '',
           '1989 - 1994 (n = 100)', '', '1995 - 1998 (n = 114)', ''])


# plt.xlabel('$log(Qw)$ $(m^{{3}}yr^{-1})$')
# plt.ylabel('$log(Qs)$ $(m^{3}yr^{-1})$')
plt.xlabel('log(water discharge) $(m^{{3}}$ $yr^{-1})$')
plt.ylabel('log(sediment discharge) $(m^{3}$ $yr^{-1})$')
plt.text(20.5, 4, 'y = ' + str(m) + 'x ' + str(b))
plt.grid()
plt.tight_layout()
fig.savefig(str(Directory)+'/Output/Rating_Curve_Grouped.' +
            Export_format, dpi=dpi)
plt.close(fig)


# %% Use rating curve to calculate annual sediment export

# Initiate variables
Qw = []
Qs = []
Qs_annual = []
Qs_max_annual = []
max_Qs = 0
#ticker = 0
Qs_sum = 0

# Execute
# for i in np.arange(92, np.size(Gauge_Hydro['day'])):
#     Qw = Gauge_Hydro['daily_discharge'][i]
#     Qw = np.log(Qw)
#     Qs = (Qw * m) + b
#     Qs = e ** Qs
#     Qs_sum += Qs
#     if Qs > max_Qs:
#         max_Qs = Qs
#     ticker += 1
#     if ticker >= 365:
#         Qs_annual = np.append(Qs_annual, Qs_sum)
#         Qs_max_annual = np.append(Qs_max_annual, max_Qs)
#         Qs_sum = 0
#         max_Qs = 0
#         ticker = 0
# mcer = Qs_annual / gauge_upstream_area

for i in np.arange(0, np.size(Hydro['daily_discharge_m3yr-1'])):
    Qw = Hydro['daily_discharge_m3yr-1'][i]
    Qw = np.log(Qw)
    Qs = (Qw * m) + b
    Qs = e ** Qs

    # Account for time
    Qs *= (1 / 365)     # dt = 1 / 365

    Qs_sum += Qs
    if Qs > max_Qs:
        max_Qs = Qs
    #ticker += 1
    if i < np.size(Hydro['daily_discharge_m3yr-1']) - 1:
        if Hydro['year'][i + 1] > Hydro['year'][i]:
            Qs_annual = np.append(Qs_annual, Qs_sum)
            Qs_max_annual = np.append(Qs_max_annual, max_Qs)
            Qs_sum = 0
            max_Qs = 0
            #ticker = 0
    if i == np.size(Hydro['daily_discharge_m3yr-1']) - 1:
        Qs_annual = np.append(Qs_annual, Qs_sum)
        Qs_max_annual = np.append(Qs_max_annual, max_Qs)
        Qs_sum = 0
        max_Qs = 0
        #ticker = 0
mcer = Qs_annual / gauge_upstream_area
mcer_mean = np.nanmean(mcer)
mcer_std = np.nanstd(mcer)
mcer_25 = np.nanpercentile(mcer, 25)
mcer_75 = np.nanpercentile(mcer, 75)


# %% Plot estimated annual sediment exported

if path.exists(str(Directory)+'/Output') == False:
    os.mkdir(str(Directory)+'/Output')
plt.ioff()
fig = plt.figure()
ax = fig.add_subplot()
#xaxis = np.arange(2000, 2021)
xaxis = np.arange(1930, 2022)
xaxis = xaxis.astype(str)
plt.plot(xaxis, Qs_annual, 'ko')
xaxis[1: -1] = ''
ax.set_xticklabels(xaxis, rotation=45)
plt.xlabel('Year')
plt.ylabel('Estimated annual sediment flux $(m^{3})$')
plt.grid()
plt.tight_layout()
fig.savefig(str(Directory) +
            '/Output/Estimated_Annual_Exported_Sediment.' + Export_format, dpi=dpi)
plt.close(fig)


# %% Plot relationship between estimated annual export and largest event of each year

# if path.exists(str(Directory)+'/Output') == False:
#     os.mkdir(str(Directory)+'/Output')
# plt.ioff()
# fig = plt.figure()
# plt.plot(Qs_max_annual, Qs_annual, 'o')
# plt.xlabel('Predicted sediment exported during largest annual flow event (m^3)')
# plt.ylabel('Predicted total exported sediment (m^3)')
# plt.grid()
# plt.tight_layout()
# fig.savefig(str(Directory) +
#             '/Output/Relationship_With_Largest_Annual_Event.' + Export_format, dpi=dpi)
# plt.close(fig)


# %% Same as above but as a percentage

# if path.exists(str(Directory)+'/Output') == False:
#     os.mkdir(str(Directory)+'/Output')
# plt.ioff()
# fig = plt.figure()
# ax = fig.add_subplot()
# #xaxis = np.arange(2000, 2021)
# xaxis = np.arange(1930, 2022)
# xaxis = xaxis.astype(str)
# yaxis = (Qs_max_annual / Qs_annual) * 100
# plt.plot(xaxis, yaxis, '-o')
# xaxis[1: -1] = ''
# ax.set_xticklabels(xaxis, rotation=45)
# plt.xlabel('Year')
# plt.ylabel('% contributed by largest annual event')
# plt.grid()
# plt.tight_layout()
# fig.savefig(str(Directory) +
#             '/Output/Relationship_With_Largest_Annual_Event_Percent.' + Export_format, dpi=dpi)
# plt.close(fig)


# %% Plot predicted catchment averaged erosion rates


# if path.exists(str(Directory)+'/Output') == False:
#     os.mkdir(str(Directory)+'/Output')
# plt.ioff()
# fig = plt.figure()
# ax = fig.add_subplot()
# #xaxis = np.arange(2000, 2021)
# xaxis = np.arange(1930, 2022)
# xaxis = xaxis.astype(str)
# plt.plot(xaxis, mcer, '-o', color='b', label='$Annual$ $predictions$')
# plt.axhline(mcer_mean, color='b', linestyle='--', label='$Mean$')
# #plt.axhspan(mcer_25, mcer_mean + mcer_std, alpha = 0.5, color = 'b')
# #plt.axhspan(SAP56_mean_recalc - SAP56_std_recalc, SAP56_mean_recalc + SAP56_std_recalc, alpha = 0.5, color = 'r', label = '$^{10}Be$ $implied$')
# #plt.axhspan(SAP60_mean - SAP60_std, SAP60_mean + SAP60_std, alpha = 0.5, color = 'm')
# #plt.axhline(Trimble_yield, color = 'g')
# #plt.axhline(Trimble_upland, color = 'r')
# xaxis[1: -1] = ''
# ax.set_xticklabels(xaxis, rotation=45)
# plt.yscale('log')
# plt.xlabel('Year')
# plt.ylabel('Mean catchment erosion rate $(m$ $yr^{-1})$', fontsize=9)
# plt.grid()
# #plt.ylim([1E-6, 2E-1])
# #plt.legend(['Gage predictions', 'Gauge mean', '$^1$$^0$Be (Chestatee)', '$^1$$^0$Be (Blue Ridge, uncorrected)', '1909 sediment yield (Chattahoochie at West Point)', 'Post-industrial upland erosion (Chattahoochie above West Point)'], loc ="upper left", fontsize = 5)
# plt.legend()
# plt.tight_layout()
# fig.savefig(str(Directory) +
#             '/Output/Predicted_Catchment_Averaged_Erosion_Rates.' + Export_format, dpi=dpi)
# plt.close(fig)





#%% Compute Values to Report

Qs_annual_mean = np.nanmean(Qs_annual)              # Mean annual sediment volume exported from watershed (m^3)
Qs_annual_min = np.nanquantile(Qs_annual, 0.25)     # Min annual sediment volume exported from watershed (m^3). Defined as 25th percentile, because std(Qs_annual) > mean(Qs_annual)
Qs_annual_max = np.nanquantile(Qs_annual, 0.75)     # Max annual sediment volume exported from watershed (m^3). Defined as 75th percentile, because std(Qs_annual) > mean(Qs_annual)
E = mcer_mean
E_min = np.nanquantile(mcer, 0.25) 
E_max = np.nanquantile(mcer, 0.75) 

# Sediment Export
if path.exists(str(Directory)+'/Output') == False:
    os.mkdir(str(Directory)+'/Output')
plt.ioff()
fig = plt.figure()
ax = fig.add_subplot()
xaxis = np.arange(1930, 2022)
xaxis = xaxis.astype(str)
plt.plot(xaxis, Qs_annual, 'ko', label='Estimated annual total')
plt.axhspan(Qs_annual_min, Qs_annual_max, alpha = 0.5, color = 'k', label='Interquartile range')
plt.axhline(Qs_annual_mean, color='r', linestyle='-', label='Mean')
xaxis[1: -1] = ''
ax.set_xticklabels(xaxis, rotation=45)
plt.yscale('log')
plt.ylim(1E2, 1E5)
plt.xlabel('Year')
plt.ylabel('Annual sediment export $(m^{3})$')
plt.grid()
plt.legend()
plt.tight_layout()
fig.savefig(str(Directory) +
            '/Output/Estimated_sediment_export.' + Export_format, dpi=dpi)
plt.close(fig)

# Mean catchment erosion rate
if path.exists(str(Directory)+'/Output') == False:
    os.mkdir(str(Directory)+'/Output')
plt.ioff()
fig = plt.figure()
ax = fig.add_subplot()
xaxis = np.arange(1930, 2022)
xaxis = xaxis.astype(str)
plt.plot(xaxis, mcer, 'ko', label='Estimated annual rate')
plt.axhspan(E_min, E_max, alpha = 0.5, color = 'k', label='Interquartile range')
plt.axhline(E, color='g', linestyle='-', label='Mean')
plt.axhline(np.nanmedian(mcer), color='orange', linestyle='-', label='Median')
xaxis[1: -1] = ''
ax.set_xticklabels(xaxis, rotation=45)
plt.yscale('log')
plt.ylim(1E-6, 1E-3)
plt.xlabel('Year')
plt.ylabel('Mean catchment erosion rate $(m^{3} yr^{-1})$')
plt.grid()
plt.legend()
plt.tight_layout()
fig.savefig(str(Directory) +
            '/Output/Estimated_catchment_averaged_erosion_rates.' + Export_format, dpi=dpi)
plt.close(fig)

# Past and Present Mean catchment erosion rate
if path.exists(str(Directory)+'/Output') == False:
    os.mkdir(str(Directory)+'/Output')
plt.ioff()
fig = plt.figure()
ax = fig.add_subplot()
xaxis = np.arange(1930, 2022)
xaxis = xaxis.astype(str)
plt.plot([], [], marker='None', linestyle='None', label =  r"$\bf{" + 'Rating' + "}$" + ' ' + r"$\bf{" + 'curve' + "}$")
plt.plot(xaxis, mcer, 'ko', label='Estimated annual rates')
plt.axhline(E, color='g', linestyle='-', label='Mean')
plt.axhline(np.nanmedian(mcer), color='orange', linestyle='-', label='Median')
plt.axhspan(E_min, E_max, alpha = 0.5, color = 'k', label='Interquartile range')
plt.plot([], [], marker='None', linestyle='None', label= r'$\bf{^{10}Be}$' )
plt.axhline(7.96E-06, color='b', linestyle='-', label='mean')
plt.axhspan(6.19E-06, 9.73E-06, alpha = 0.5, color = 'c', label='+/- 1 sigma')
plt.plot([], [], marker='None', linestyle='None', label=' ' )
xaxis[1: -1] = ''
ax.set_xticklabels(xaxis, rotation=45)
plt.yscale('log')
plt.ylim(1E-6, 1E-3)
plt.xlabel('Year')
plt.ylabel('Mean catchment erosion rate $(m^{3} yr^{-1})$')
plt.grid()
plt.legend(loc=2, ncol=2, prop={'size': 6})
plt.tight_layout()
fig.savefig(str(Directory) +
            '/Output/Both_Catchment_Averaged_Erosion_Rates.' + Export_format, dpi=dpi)
plt.close(fig)

# E1 vs E1rc
if path.exists(str(Directory)+'/Output') == False:
    os.mkdir(str(Directory)+'/Output')
plt.ioff()
fig, ax = plt.subplots()
plt.boxplot(mcer[~np.isnan(mcer)], showmeans=True, meanline=True)
plt.errorbar(1.3, 2.94694647344916E-05, 6.55288348995580E-06, capsize=2, ecolor='red')
# plt.errorbar(1.6, 7.96E-06, 1.77E-06, capsize=2, ecolor='b')
plt.ylabel('Mean catchment erosion rate $(m^{3} yr^{-1})$')
fig.set_figwidth(2)
plt.tight_layout()
ax = plt.gca()
ax.xaxis.set_tick_params(labelbottom=False)
ax.set_xticks([])
plt.ylim(1E-6, 1E-3)
plt.yscale('log')
plt.grid(which='both', linestyle='--', linewidth=0.5)
# plt.show()
fig.savefig(str(Directory) + '/Output/E1_E1rc_comparison.' + Export_format, dpi=dpi)
plt.close(fig)

# E1 vs E1rc with 10Be
if path.exists(str(Directory)+'/Output') == False:
    os.mkdir(str(Directory)+'/Output')
plt.ioff()
fig, ax = plt.subplots()
plt.boxplot(mcer[~np.isnan(mcer)], showmeans=True, meanline=True)
plt.errorbar(1.3, 2.94694647344916E-05, 6.55288348995580E-06, capsize=2, ecolor='red')
plt.errorbar(1.6, 7.96E-06, 1.77E-06, capsize=2, ecolor='b')
plt.errorbar(3, 7.96E-06, 1.77E-06, alpha=0)
plt.ylabel('Mean catchment erosion rate $(m^{3} yr^{-1})$')
fig.set_figwidth(3)
plt.tight_layout()
ax = plt.gca()
ax.xaxis.set_tick_params(labelbottom=False)
ax.set_xticks([])
plt.ylim(1E-6, 1E-3)
plt.yscale('log')
plt.grid(which='both', linestyle='--', linewidth=0.5)
# plt.show()
fig.savefig(str(Directory) + '/Output/E1_E1rc_10Be_comparison.' + Export_format, dpi=dpi)
plt.close(fig)





#%% Combined

# Mean catchment erosion rate
if path.exists(str(Directory)+'/Output') == False:
    os.mkdir(str(Directory)+'/Output')
plt.ioff()
#
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
#
xaxis = np.arange(1930, 2022)
xaxis = xaxis.astype(str)
ax1.plot(xaxis, mcer, 'ko', label='Estimated annual rate')
ax1.axhspan(E_min, E_max, alpha = 0.25, color = 'k', label='Interquartile range')
ax1.axhline(E, color='g', linestyle='-', label='Mean')
ax1.axhline(np.nanmedian(mcer), color='orange', linestyle='-', label='Median')
xaxis[1: -1] = ''
ax1.set_xticklabels(xaxis, rotation=45)
plt.yscale('log')
plt.ylim(1E-6, 1E-3)
ax1.set_xlabel('Year')
ax1.set_ylabel('Mean catchment erosion rate $(m^{3} yr^{-1})$')
# ax1.grid()
ax1.grid(which='both', axis='y', linestyle='--', linewidth=0.5)
ax1.legend(prop={'size': 6})
plt.tight_layout()
#
bp1 = ax2.boxplot(mcer[~np.isnan(mcer)], showmeans=True, meanline=True, patch_artist=True, boxprops=dict(facecolor="k", alpha = 0.25))
eb1 = ax2.errorbar(1.3, 2.94694647344916E-05, 6.55288348995580E-06, capsize=2, color='red', ecolor='red', label='Test2')
eb2 = ax2.errorbar(1.6, 7.96E-06, 1.77E-06, capsize=2, color='blue', ecolor='b', label='Test3')
# ax2.errorbar(3, 7.96E-06, 1.77E-06, alpha=0)
# plt.ylabel('Mean catchment erosion rate $(m^{3} yr^{-1})$')
# fig.set_figwidth(3)
# plt.tight_layout()
# ax = plt.gca()
ax2.xaxis.set_tick_params(labelbottom=False)
# ax.set_xticks([])
# ax2.set_xticklabels(xaxis, rotation=45)
# plt.ylim(1E-6, 1E-3)
# plt.yscale('log')
ax2.grid(which='both', axis='y', linestyle='--', linewidth=0.5)
#ax2.legend()
ax2.legend([bp1["boxes"][0], eb1, eb2], ['Rating curve data', '$E_{1}$', '$^{10}Be$'], loc='upper right', prop={'size': 6})
#
fig.savefig(str(Directory) +
            '/Output/Combined.' + Export_format, dpi=dpi)
plt.close(fig)




# Rating_Curve







