#!/usr/bin/env python
# coding: utf-8

# # TerrainSandbox 
# 
# **Author**: Chris Sheehan

#%% Block 1: The Grid

import inspect
import os
script_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
project_directory = script_directory.replace('\k1u_empirical', '')

# ENTER VARIABLES ############################################################
##############################################################################

m_sp = 0.503                         
n_sp = 1.224                         
K = 6.16968E-07
D = (K * 800) 

# DEM_path = 'C:/Users/sheehacz/Dropbox/BC_Landlab/SPACE_2022/New_Models/Input/Chestatee.asc'        # '/Users/Chris/Dropbox/BC_Landlab/SPACE_2022/Chestatee/D_Calibration_New/DEMs/A.asc'
DEM_path = project_directory + '\DEM\Chestatee.asc'
dxy = 26.673674998310        
no_data_value = -99999   
pour_point_node = 327       # pre-selected. Ideal beacause (A). Running  FlowAccumulator on raw DEM identifies this node, and (B). It is a boundary node. Chestatee = 327. A = 3959.

#lc_path = 'C:/Users/sheehacz/Dropbox/BC_Landlab/SPACE_2022/New_Models/Input/NCLD_2019_Align.asc'
#lc_scheme_path = 'C:/Users/sheehacz/Dropbox/BC_Landlab/SPACE_2022/New_Models/Input/Erosion_Rate_Scheme.csv'

SAP56_mean_recalc = 7.96 * 1E-6          # m / yr^-1
SAP56_mean_recalc_std = 1.77 * 1E-6

Total_area = 608568641.6459954        # m^2, calculated using flow router

Directory = script_directory + '/'
Export_format = 'png'
dpi = 600

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
from landlab.components import PerronNLDiffuse, StreamPowerEroder, Space, DepressionFinderAndRouter, LinearDiffuser, TaylorNonLinearDiffuser,  FlowAccumulator, ChannelProfiler, SteepnessFinder, ChiFinder, Lithology, LithoLayers, NormalFault
from landlab.io.esri_ascii import write_esri_ascii
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import colors as colors
from landlab.plot.graph import plot_graph
from landlab.io.esri_ascii import read_esri_ascii 
import EET
import matplotlib.cbook as cbook

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
mg.add_zeros('node', 'soil__flux')

# Handle Grid Boundaries
print('Handling grid boundaries...')
mg.set_status_at_node_on_edges(right=4, top=4, left=4, bottom=4)
mg.status_at_node[pour_point_node] = mg.BC_NODE_IS_FIXED_VALUE

# Enforce domain as a single watershed
print('Enforcing domain as a single watershed...')
number_of_watersheds = 1

# Record pour point node elevation
print('Identifying pour point...')
pour_point_node_elevation = zr[pour_point_node]

# Set uplift rate
print('Setting uplift rate...')
U = np.ones(mg.number_of_nodes) * SAP56_mean_recalc

# Initialize FlowAccumulator
print('Initializing FlowAccumulator...')
frr = FlowAccumulator(mg, flow_director='D8')
frr.run_one_step()

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

# Plot initial channel network
print('Plotting initial channel network...')
prf = ChannelProfiler(mg, number_of_watersheds=number_of_watersheds, main_channel_only=False, minimum_channel_threshold=10000000)
prf.run_one_step()
if path.exists(str(Directory)+'/TerrainSandbox/Initial_Condition') == False:
    os.mkdir(str(Directory)+'/TerrainSandbox/Initial_Condition')
plt.ioff()
fig = plt.figure()
prf.plot_profiles_in_map_view()
plt.tight_layout()
fig.savefig(str(Directory)+'/TerrainSandbox/Initial_Condition/Channel_Network.' + Export_format,  format=Export_format,  dpi=dpi)
plt.close(fig)

# # Import Land Cover
# print('Importing NLCD data...')
# lc = read_esri_ascii(lc_path)
# lc = lc[1]

# # Import land cover scheme
# print('Importing land cover scheme...')
# lc_scheme = pd.read_csv(lc_scheme_path)

# Print space
print(' ')





#%% BLOCK 5: RESET THE MODEL TIME AND THE PLOTTING TIMERS TO 0

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------

# Reset
print('Reseting time parameters...')
total_time = 0 
Plot_Ticker = 0
Export_DEM_Ticker = 0
timestep_integer = 0

# Print space
print(' ')





#%% BLOCK 6: TIME PARAMETERS

# ENTER VARIABLES ############################################################
##############################################################################

dt = 0.1          
tmax = 100 

##############################################################################
# ENTER VARIABLES ############################################################

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------

# Set t
print('Setting t...')
t = np.arange(0, tmax, dt) 

# Print space
print(' ')





#%% BLOCK 7: PLOTTING OPTIONS

# ENTER VARIABLES ############################################################
##############################################################################

# Intervals
Plot_interval = 10
Export_DEM_Interval = 1

# Drainage information
min_drainage_area = (dxy**2) * 10000
main_channel_only = False

# Toggle maps
DEM_Image = True
Channel_Map = False
Ksn_Map = False
Chi_Map = False
Erosion_Rate_Map = True
DZ_DT_Map = True
Net_Erosion_Map = True
Net_DZ_Map = True

# Toggle profiles
Slope_Area = False
Channel_Profile = False
Ksn_Profile = False
Chi_Profile = False
Erosion_Rate_Profile = False
DZ_DT_Profile = False

# Toggle timeseries
Mean_Basin_Erosion_Rate = True

# Toggle DEM export
Export_DEM = True


##############################################################################
# ENTER VARIABLES ############################################################

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------

Plot_Ticker = 0
Export_DEM_Ticker = 0

# Notify
print('Setting plot parameters')

# Print space
print(' ')





#%% BLOCK 8: TIME LOOP

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------

# Create directory
print('Creating directory...')
if path.exists(str(Directory)+'/TerrainSandbox') == False:
    os.mkdir(str(Directory)+'/TerrainSandbox')

# Initialize linear diffuser
print('Initializing LinearDiffuser...') 
dfn = LinearDiffuser(mg, 
                      linear_diffusivity = D, 
                      method = 'simple',
                      deposit = False
                      )        

# Initialize StreamPowerEroder 
print('Initializing StreamPowerEroder...')  
spr = StreamPowerEroder(mg, 
                        K_sp = K, 
                        m_sp = m_sp, 
                        n_sp = n_sp, 
                        erode_flooded_nodes = True
                        )

# Initialize DepressionFinderAndRouter
print('Initializing DepressionFinderAndRouter...') 
dfr = DepressionFinderAndRouter(mg)

# Initialize ErosionElevationTracker
print('Initializing ErosionElevationTracker...') 
eet = EET.ErosionElevationTracker(mg, bedrock__and__soil = False)

# Initialize appendable arrays
print('Initializing appendable arrays') 
mean_basin_erosion_rate = []
sum_basin_erosion = []
timesteps = []
times = []
#
std_basin_erosion_rate = []
stats_list = list()

# Initialize ChannelProfiler
print('Initializing ChannelProfiler') 
if main_channel_only == True:
    prf = ChannelProfiler(mg, number_of_watersheds=number_of_watersheds, main_channel_only=True, minimum_channel_threshold=min_drainage_area)
if main_channel_only == False:
    prf = ChannelProfiler(mg, number_of_watersheds=number_of_watersheds, main_channel_only=False, minimum_channel_threshold=min_drainage_area)

# Print space
print(' ')





#%% Calculate

E = 2.94694647344916E-05
num_nodes = np.size(zr) - np.size(no_data_nodes)

Ai_m = mg.at_node['drainage_area'] ** m_sp
Ai_m[no_data_nodes] = nan
Si_n = mg.at_node['topographic__steepest_slope'] ** n_sp
Si_n[no_data_nodes] = nan

mult = Ai_m * Si_n
add = np.nansum(mult)
divd = add / num_nodes

K_empir = E / divd
K_empir_min = (2.29165812445355E-05) / divd
K_empir_max = (3.60223482244474E-05) / divd






















  