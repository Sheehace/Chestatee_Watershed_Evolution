#!/usr/bin/env python
# coding: utf-8

# # TerrainSandbox 
# 
# **Author**: Chris Sheehan

#%% Block 1: The Grid

import inspect
import os
script_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
project_directory = script_directory.replace('\k1_create', '')

# ENTER VARIABLES ############################################################
##############################################################################

# DEM_path = 'C:/Users/sheehacz/Dropbox/BC_Landlab/SPACE_2022/New_Models/Input/Chestatee.asc'      # '/Users/Chris/Dropbox/BC_Landlab/SPACE_2022/Chestatee/DEMs/Chestatee.asc'
DEM_path = project_directory + '\DEM\Chestatee.asc'
dxy = 26.673674998310        
no_data_value = -99999   
pour_point_node = 327       # pre-selected. Ideal beacause (A). Running  FlowAccumulator on raw DEM identifies this node, and (B). It is a boundary node. Chestatee = 327. A = 3959.

# lc_path = 'C:/Users/sheehacz/Manuscript_Codes/K1_Create/NCLD_2019_Align.asc'
# lc_scheme_path = 'C:/Users/sheehacz/Manuscript_Codes/K1_Create/Erosion_Rate_Scheme.csv'
lc_path = script_directory + '/NCLD_2019_Align.asc'
lc_scheme_path = script_directory + '/Erosion_Rate_Scheme.csv'

SAP56_mean_recalc = 7.96 * 1E-6          # m / yr^-1
SAP56_std_recalc = 0.63 * 1E-6        # m / yr^-1
SAP56_mean = 9.5 * 1E-6
SAP56_std = 0.75 * 1E-6
SAP60_mean = 17.16 * 1E-6
SAP60_std = 1.32 * 1E-6

mcer_mean = 3.331058457066812e-05 # 0.01956698905849651  # 0.018050098405579439    # m / yr

Total_area = 608568641.6459954        # m^2, calculated using flow router
Kb_area = 119163054.826772                  # Over_590 = 119163054.826772

Directory = script_directory + '/'
Export_format = 'pdf'
dpi = 600

Burn_min_drainage_area = (dxy**2) * 100       # (dxy**2) * 10000
Burn_depth_min = 0.1
Burn_depth_max = 10

##############################################################################
# ENTER VARIABLES ############################################################

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------

# Print space
print(' ')

# Import libraries
print('Importing libraries...')
import numpy as np
from numpy import nan, isnan
import pandas as pd
import math
import os
import os.path
from os import path
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
warnings.filterwarnings("ignore")
from landlab import RasterModelGrid, imshow_grid
from landlab.plot.video_out import VideoPlotter
from landlab.components import PerronNLDiffuse, StreamPowerEroder, Space, SpaceLargeScaleEroder, DepressionFinderAndRouter, LinearDiffuser, TaylorNonLinearDiffuser,  FlowAccumulator, ChannelProfiler, SteepnessFinder, ChiFinder, Lithology, LithoLayers, NormalFault
from landlab.io.esri_ascii import write_esri_ascii
from matplotlib import pyplot as plt
from matplotlib import cm
from landlab.plot.graph import plot_graph
from landlab.io.esri_ascii import read_esri_ascii 

# Set random seed
print('Setting random seed...')
np.random.seed(0) 

# Create directories
print('Creating directory...')
if path.exists(str(Directory)+'/TerrainSandbox') == False:
    os.mkdir(str(Directory)+'/TerrainSandbox')
if path.exists(str(Directory)+'/TerrainSandbox/Initial_Condition') == False:
    os.mkdir(str(Directory)+'/TerrainSandbox/Initial_Condition')
if path.exists(str(Directory)+'/TerrainSandbox/Analysis') == False:
    os.mkdir(str(Directory)+'/TerrainSandbox/Analysis')
if path.exists(str(Directory)+'/TerrainSandbox/CSV') == False:
    os.mkdir(str(Directory)+'/TerrainSandbox/CSV')

# Import DEM
print('Importing DEM...')
(mg, zr) = read_esri_ascii(DEM_path, name='topographic__elevation')

# Create copy of DEM
print('Copying DEM for initial topography comparisons...')
(mg0, zr0) = read_esri_ascii(DEM_path, name='topographic__elevation')

# Handle model DEM dimensions and non-value nodes
print('Handling DEM dimensions and non-value nodes...')
mg.set_nodata_nodes_to_closed(zr, no_data_value)
mg0.set_nodata_nodes_to_closed(zr0, no_data_value)
no_data_nodes = np.where(mg.at_node['topographic__elevation'] == no_data_value)
no_data_nodes = no_data_nodes[0]
number_of_rows = mg.number_of_cell_rows
number_of_columns = mg.number_of_cell_columns
nodes = np.arange(0, np.size(zr))

# Create node keys
print('Creating node keys...')
erosion_rate = mg.add_zeros('node', 'erosion_rate')
#erode_depo = mg.add_zeros('node', 'erode_depo')
previous_zr = mg.add_zeros('node', 'previous_topographic__elevation')
dz_dt = mg.add_zeros('node', 'elevational_change')
mg.add_zeros('node', 'soil__flux')

# Handle Grid Boundaries
mg.set_status_at_node_on_edges(right=4, top=4, left=4, bottom=4)
mg.status_at_node[pour_point_node] = mg.BC_NODE_IS_FIXED_VALUE

# Record pour point node elevation
pour_point_node_elevation = zr[pour_point_node]

# Burn stream network #########################################################

# # Display status
# print('Burning DEM...')

# Initialize FlowAccumulator
frr = FlowAccumulator(mg, flow_director='D8')
frr.run_one_step()

# # Initialize ChannelProfiler
# prf = ChannelProfiler(mg, number_of_watersheds=1, main_channel_only=False, minimum_channel_threshold=Burn_min_drainage_area)
# prf.run_one_step()

# # Find channel nodes and their drainage areas
# ids = []
# for i, outlet_id in enumerate(prf.data_structure):
#     for j, segment_id in enumerate(prf.data_structure[outlet_id]):
#         segment = prf.data_structure[outlet_id][segment_id]
#         profile_ids = segment["ids"]
#         ids = np.append(ids, profile_ids)        
# ids = ids.astype(int)
# ids_DA = mg.at_node['drainage_area'][ids]

# # Calculate burn
# Burn_depth = ((ids_DA / np.max(ids_DA)) * Burn_depth_max) + Burn_depth_min

# # Flatten terrain 
# zr[np.where(mg.at_node['topographic__elevation'] > no_data_value)] = pour_point_node_elevation + Burn_depth_max

# # Execute burn
# mg.at_node['topographic__elevation'][ids] -= Burn_depth

# # Add topographic noise to terrain outside of channels
# flats = np.where(zr == np.max(zr))
# flats = flats[0]
# mg_noise = np.random.rand(np.size(flats))/1000. 
# zr[flats] += mg_noise

# # Method B
# zr = ( ((zr - pour_point_node_elevation) / (np.max(zr) - pour_point_node_elevation)) * 100 ) + pour_point_node_elevation
# zr[no_data_nodes] = no_data_value
# mg.at_node['topographic__elevation'] = zr


# Burn stream network #########################################################

# Plot initial model terrain
print('Plotting initial model terrain...')
if path.exists(str(Directory)+'/TerrainSandbox/Initial_Condition') == False:
    os.mkdir(str(Directory)+'/TerrainSandbox/Initial_Condition')
plt.ioff()
fig = plt.figure()         
imshow_grid(mg, 'topographic__elevation', grid_units=('m', 'm'), var_name="Elevation (m)", cmap='terrain', allow_colorbar=True)
plt.tight_layout()
fig.savefig(str(Directory)+'/TerrainSandbox/Initial_Condition/DEM_Image.' + Export_format, dpi = dpi)
plt.close(fig)

# # Plot initial channel network
# print('Plotting initial channel network...')
# prf.run_one_step()
# if path.exists(str(Directory)+'/TerrainSandbox/Initial_Condition') == False:
#     os.mkdir(str(Directory)+'/TerrainSandbox/Initial_Condition')
# plt.ioff()
# fig = plt.figure()
# prf.plot_profiles_in_map_view()
# plt.tight_layout()
# fig.savefig(str(Directory)+'/TerrainSandbox/Initial_Condition/Channel_Map.' + Export_format,  format=Export_format,  dpi = dpi)
# plt.close(fig)

# # Import K_zones
# print('Importing K_zone...')
# (mgk, K_zone) = read_esri_ascii(K_zone_path, name='temp')
# K_zone[np.where(K_zone == np.min(K_zone))] = nan
# mg.add_ones('node', 'K_zone')
# mg.at_node['K_zone'] *= K_zone
# del(mgk, K_zone)

# Rate equations in notebook ##################################################

# Set parameters
Ka_area = Total_area - Kb_area
f_a = Ka_area / Total_area
f_b = Kb_area / Total_area
rate_b = SAP60_mean
rate_total = SAP56_mean

# Calculate rate_a
rate_a = ( rate_total - (rate_b * f_b) ) / f_a

# Calculate fraction Kb in Chestatee Basin
Kb_Ka_estimate = rate_b / rate_a

# Calculate Ka erosion rate implied by 

# Rate equations in notebook ##################################################

# Import Land Cover
lc = read_esri_ascii(lc_path)
lc = lc[1]

# Import land cover scheme
lc_scheme = pd.read_csv(lc_scheme_path)

# if path.exists('C:/Users/sheehacz/Manuscript_Codes/K1_Create/Output') == False:
#     os.mkdir('C:/Users/sheehacz/Manuscript_Codes/K1_Create/Output')
if path.exists(script_directory + '/Output') == False:
    os.mkdir(script_directory + '/Output')


# Create K
print("Setting K")
Ka_base = 8.29124436168032E-07
K = np.ones(mg.number_of_nodes) * nan
for i in np.arange(0, np.size(lc)):
    if zr[i] != no_data_value:
        value = lc[i]
        r = np.where(lc_scheme['nlcd_value'] == value)
        r = r[0]
        K[i] = lc_scheme['Normalized Prefered C Factor (Ci)'][r] * Ka_base
np.savetxt(script_directory + '/Output/K1.csv', K, delimiter = ",")

# Script below modified from K0_Empirical_pc.py
num_nodes = np.size(zr) - np.size(no_data_nodes)

m_sp = 0.503                         
n_sp = 1.224                          

Ai_m = mg.at_node['drainage_area'] ** m_sp
Ai_m[no_data_nodes] = nan
Si_n = mg.at_node['topographic__steepest_slope'] ** n_sp
Si_n[no_data_nodes] = nan

mult = Ai_m * Si_n * K
E1_mean = np.nanmean(mult)

Kmin = K / Ka_base * 6.44758826618105E-07
mult = Ai_m * Si_n * Kmin
E1_min = np.nanmean(mult)

Kmax = K / Ka_base * 1.01349004571795E-06
mult = Ai_m * Si_n * Kmax
E1_max = np.nanmean(mult)

# Export K1 Map
print('Plotting initial model terrain...')
if path.exists(script_directory + '/Output') == False:
    os.mkdir(script_directory + '/Output')
plt.ioff()
fig = plt.figure()         
imshow_grid(mg, K, grid_units=('m', 'm'), var_name="Erodibility", cmap='plasma', allow_colorbar=True)
plt.tight_layout()
fig.savefig(script_directory + '/Output/K1.' + Export_format, dpi = dpi)
plt.close(fig)

# Focus area 1
# axes = plt.gca()
# axes.get_ylim()
#
plt.ioff()
fig = plt.figure()
imshow_grid(mg, K, cmap="plasma")
plt.xlim([775630, 781890]) # 6260
plt.ylim([3820650, 3826330]) # 5680
plt.tight_layout()
plt.show()
fig.savefig(script_directory + '/Output/Focus_1' + '.' + Export_format, format=Export_format, dpi=dpi)
# fig.savefig(script_directory + '/Output/Focus_1' + '.png', format='png', dpi=dpi)
plt.close(fig)

# Focus area 2
# axes = plt.gca()
# axes.get_ylim()
#
plt.ioff()
fig = plt.figure()
imshow_grid(mg, K, cmap="plasma")
plt.xlim([775050, 783490])
plt.ylim([3840900, 3848560]) # 7660
plt.tight_layout()
plt.show()
fig.savefig(script_directory + '/Output/Focus_2' + '.' + Export_format, format=Export_format, dpi=dpi)
# fig.savefig(script_directory + '/Output/Focus_2' + '.png', format='png', dpi=dpi)
plt.close(fig)

# Focus area 3
# axes = plt.gca()
# axes.get_ylim()
#
plt.ioff()
fig = plt.figure()
imshow_grid(mg, K, cmap="plasma")
plt.xlim([788790, 797220])      # 8880
plt.ylim([3823370, 3831020])    # 7650
plt.tight_layout()
plt.show()
fig.savefig(script_directory + '/Output/Focus_3' + '.' + Export_format, format=Export_format, dpi=dpi)
# fig.savefig(script_directory + '/Output/Focus_3' + '.png', format='png', dpi=dpi)
plt.close(fig)
