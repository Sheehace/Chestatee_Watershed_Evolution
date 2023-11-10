#!/usr/bin/env python
# coding: utf-8

# # TerrainSandbox 
# 
# **Author**: Chris Sheehan

#%% Block 1: The Grid

import inspect
import os
script_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
project_directory = script_directory.replace('\d', '')

# ENTER VARIABLES ############################################################
##############################################################################

#DEM_path = 'C:/Users/sheehacz/Dropbox/BC_Landlab/SPACE_2022/New_Models/Input/Chestatee.asc'      # '/Users/Chris/Dropbox/BC_Landlab/SPACE_2022/Chestatee/DEMs/Chestatee.asc'
DEM_path = project_directory + '\DEM\Chestatee.asc'

dxy = 26.673674998310        
no_data_value = -99999   
pour_point_node = 327       # pre-selected. Ideal beacause (A). Running  FlowAccumulator on raw DEM identifies this node, and (B). It is a boundary node. Chestatee = 327. A = 3959, B = 98, C = 7738. 

# K_zone_path = '/Users/Chris/Dropbox/BC_Landlab/SPACE_2022/Chestatee/Land_Cover/K_zone/Over_590.asc'
# Kb_Ka = 2

SAP56_mean_recalc = 7.96 * 1E-6          # m / yr^-1
# SAP56_std_recalc = 0.63 * 1E-6        # m / yr^-1
# SAP56_mean = 9.5 * 1E-6
# SAP56_std = 0.75 * 1E-6
# SAP60_mean = 17.16 * 1E-6
# SAP60_std = 1.32 * 1E-6

# Total_area = 608568641.6459954        # m^2, calculated using flow router
# Kb_area = 119163054.826772                  # Over_590 = 119163054.826772

Directory = script_directory + '/'
Export_format = 'png'
dpi = 600

# Burn_min_drainage_area = (dxy**2) * 100       # (dxy**2) * 10000
# Burn_depth_min = 0.1
# Burn_depth_max = 10

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

# # Initialize FlowAccumulator
# frr = FlowAccumulator(mg, flow_director='D8')
# frr.run_one_step()

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

# # Set parameters
# Ka_area = Total_area - Kb_area
# f_a = Ka_area / Total_area
# f_b = Kb_area / Total_area
# rate_b = SAP60_mean
# rate_total = SAP56_mean

# Calculate rate_a
# rate_a = ( rate_total - (rate_b * f_b) ) / f_a

# Calculate fraction Kb in Chestatee Basin
# Kb_Ka_estimate = rate_b / rate_a

# Calculate Ka erosion rate implied by 

# Rate equations in notebook ##################################################





#%% BLOCK 3: ENTER LITHOLOGIC PARAMETERS 

# ENTER VARIABLES ############################################################
##############################################################################

m_sp = 0.503                         # 0.491149684
n_sp = 1.224                          # 1.193861385

Ka = 5E-7
D = (Ka * 375) 

##############################################################################
# ENTER VARIABLES ############################################################

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------

# Print space
print(' ')

# Set K
print("Setting K")
K = np.ones(mg.number_of_nodes) * Ka
#K[np.where(mg.at_node['K_zone'] == 1)] = Ka / Kb_Ka

# Set D
# print("Setting D")
# D = np.ones(mg.number_of_nodes) * D

# Set Lith_dict
Lith_dict = {"Ksp_1": Ka}
Min_Ksp = min(Lith_dict.values())    
    
# # Plot Spatial_K
# print('Plotting spatial K...')
# if path.exists(str(Directory)+'/TerrainSandbox/Initial_Condition') == False:
#     os.mkdir(str(Directory)+'/TerrainSandbox/Initial_Condition')
# plt.ioff()
# fig = plt.figure()         
# imshow_grid(mg, K, grid_units=('m', 'm'), var_name='K', cmap='viridis', allow_colorbar=True)
# plt.tight_layout()
# fig.savefig(str(Directory)+'/TerrainSandbox/Initial_Condition/Spatial_K.' + Export_format, dpi = dpi)
# plt.close(fig)    
    




#%% BLOCK 4: TECTONIC PARAMETERS

# ENTER VARIABLES ############################################################
##############################################################################

U_1 = SAP56_mean_recalc

##############################################################################
# ENTER VARIABLES ############################################################

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------

# Print space
print(' ')

# Set U
print("Setting U")
U = np.ones(mg.number_of_nodes) * U_1

# Set U_dict
U_dict = {"U_1": U_1}
Max_U = max(U_dict.values())





#%% BLOCK 5: RESET THE MODEL TIME AND THE PLOTTING TIMERS TO 0

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------

# Print space
print(' ')

# Reset
print('Reseting time parameters...')
total_time = 0 
Plot_Ticker = 0
Export_DEM_Ticker = 0
U_Ticker = 0
Fault_Ticker = 0               





#%% BLOCK 6: TIME PARAMETERS

# ENTER VARIABLES ############################################################
##############################################################################

dt = 1E6          
tmax = 1E9 

#D = (Ka * 200) * (1E5 / dt)

##############################################################################
# ENTER VARIABLES ############################################################

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------

# Print space
print(' ')

# Set t
print('Setting t...')
t = np.arange(0, tmax, dt) 





#%% BLOCK 7: PLOTTING OPTIONS

# ENTER VARIABLES ############################################################
##############################################################################

Plot_interval = 1E8
Export_DEM_Interval = 1E9
number_of_watersheds = 1
min_drainage_area = 700 #(dxy**2) 
main_channel_only = False

DEM_Image = True
Slope_Area = False
Channel_Profile = False
Channel_Map = False
Ksn_Profile = False
Ksn_Map = False
Chi_Profile = False
Chi_Map = False
Erosion_Rate_Profile = False
Erosion_Rate_Map = False
DZ_DT_Profile = False
DZ_DT_Map = True
Ksp_Profile = False
Ksp_Map = False

Mean_Basin_Erosion_Rate = False
Zr_Diff = False

Export_DEM = True

Terrain_3D = False

# CHOOSE FROM ONE OF THESE THREE COLOR SCALING OPTIONS. IF "Color_Scaling_Prescribed" 
# IS SELECTED, ENTER A VALUE FOR "Max_color_scale"
Color_Scaling_Automatic = True
Color_Scaling_Updated = False
Color_Scaling_Prescribed = False
Max_color_scale = 10 

# CHOOSE FROM ONE OF THESE TWO ELEVATION EXAGGERATION OPTIONS. IF "Elevation_Exaggeration_Prescribed" 
# IS SELECTED, ENTER A CALUE FOR "Exaggeration_factor"
Elevation_Exaggeration_Automatic = True
Elevation_Exaggeration_Prescribed = False
Exaggeration_factor = 10

##############################################################################
# ENTER VARIABLES ############################################################

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------

Plot_Ticker = 0
Export_DEM_Ticker = 0





#%% BLOCK 8: TIME LOOP

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------

# Print space
print(' ')

# Create directory
print('Creating directory...')
if path.exists(str(Directory)+'/TerrainSandbox') == False:
    os.mkdir(str(Directory)+'/TerrainSandbox')

# Calculate ceiling
print('Calculating ceiling...')
ceiling = (Max_U / Min_Ksp) * 10

# Initialize FlowAccumulator
print('Initializing FlowAccumulator...')
frr = FlowAccumulator(mg, flow_director='D8')
frr.run_one_step()

# Initialize TaylorNonLinearDiffuser
print('Initializing TaylorNonLinearDiffuser...')
nld = TaylorNonLinearDiffuser(grid = mg,
                              linear_diffusivity = D,
                              slope_crit = 0.75315785,          # 99th percentile, topotoolbox
                              nterms = 2,
                              dynamic_dt = True,
                              if_unstable='warn',
                              courant_factor = 0.2
                              )

# nld = PerronNLDiffuse(grid = mg, 
#                       nonlinear_diffusivity = D, 
#                       S_crit = 0.75315785,                         # 99th percentile, topotoolbox
#                       )

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
dfr = DepressionFinderAndRouter(mg)
 
# Initialize previous_zr
print('Initializing previous_zr...')  
previous_zr[mg.core_nodes] = zr[mg.core_nodes]

# Initialize mean_basin_erosion_rate
mean_basin_erosion_rate = []

# Initialize comparision arrays
Diff = {'total_zr_diff_mean' : [],
        'total_zr_diff_std' : [],
        'A_zr_diff_mean' : [],
        'A_zr_diff_std' : [],
        'B_zr_diff_mean' : [],
        'B_zr_diff_std' : []
        }

# Print space
print(' ')

# Initialize previous_zr
print('Starting time loop...')  
print(' ')
for ti in t:
      
    # Uplift topograpghy
    zr[mg.core_nodes] += U[mg.core_nodes]*dt    
    
    # Run one steps                       
    frr.run_one_step()  
    dfr.map_depressions()                                  
    spr.run_one_step(dt)
    dfn.run_one_step(dt)
    #nld.run_one_step(dt)
    
    # Calculate dz_dt and erosion rate
    dz_dt[mg.core_nodes] = (zr[mg.core_nodes] - previous_zr[mg.core_nodes]) / dt
    erosion_rate[mg.core_nodes] = (previous_zr[mg.core_nodes] - (zr[mg.core_nodes] - (U[mg.core_nodes]*dt))) / dt
    
    # Record current topograpghy for comparrison before and after timestep
    previous_zr[mg.core_nodes] = zr[mg.core_nodes]    
    
    # Update mean_basin_erosion_rate
    mean_basin_erosion_rate = np.append(mean_basin_erosion_rate, np.mean(erosion_rate))
    
    # # Calculate topographic differences
    # zr_diff = mg.at_node['topographic__elevation'] - mg0.at_node['topographic__elevation'] 
    # total_zr_diff_mean = np.mean(zr_diff[mg.core_nodes])
    # total_zr_diff_std = np.std(zr_diff[mg.core_nodes])
    # A_zr_diff_mean = np.mean( zr_diff[np.where(mg.at_node['K_zone'] == 0)] )
    # A_zr_diff_std = np.std( zr_diff[np.where(mg.at_node['K_zone'] == 0)] )
    # B_zr_diff_mean = np.mean( zr_diff[np.where(mg.at_node['K_zone'] == 1)] )
    # B_zr_diff_std = np.std( zr_diff[np.where(mg.at_node['K_zone'] == 1)] )
    
    # # Append difference arrays
    # Diff['total_zr_diff_mean'] = np.append(Diff['total_zr_diff_mean'], total_zr_diff_mean)
    # Diff['total_zr_diff_std'] = np.append(Diff['total_zr_diff_std'], total_zr_diff_std)
    # Diff['A_zr_diff_mean'] = np.append(Diff['A_zr_diff_mean'], A_zr_diff_mean)
    # Diff['A_zr_diff_std'] = np.append(Diff['A_zr_diff_std'], A_zr_diff_std)
    # Diff['B_zr_diff_mean'] = np.append(Diff['B_zr_diff_mean'], B_zr_diff_mean)
    # Diff['B_zr_diff_std'] = np.append(Diff['B_zr_diff_std'], B_zr_diff_std)
    
    # Update time
    total_time += dt                                    
    print(total_time, ' ', np.max(mg.at_node['drainage_area']))
    Plot_Ticker += dt
    Export_DEM_Ticker += dt
    
    # Plots
    if Plot_Ticker == Plot_interval:
        print('Exporting figures... Please be patient!')
        
        # Initialize ChannelProfiler
        if main_channel_only == True:
            prf = ChannelProfiler(mg, number_of_watersheds=number_of_watersheds, main_channel_only=True, minimum_channel_threshold=min_drainage_area)
        if main_channel_only == False:
            prf = ChannelProfiler(mg, number_of_watersheds=number_of_watersheds, main_channel_only=False, minimum_channel_threshold=min_drainage_area)
        
        # Initialize SteepnessFinder and ChiFinder
        if total_time == Plot_Ticker:
            sf = SteepnessFinder(mg, reference_concavity=m_sp/n_sp, min_drainage_area=min_drainage_area) # NEED TO FIX! Currently breaks if you don't export images during first run of Cell 8 but then export them during later runs of Cell 8
            cf_check = 'cf' in locals()
            if cf_check == False:
                cf = ChiFinder(mg, reference_concavity=m_sp/n_sp, min_drainage_area=min_drainage_area)
        
        # DEM_Image
        if DEM_Image == True:
            if path.exists(str(Directory)+'/TerrainSandbox/DEM_Image') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/DEM_Image')
            plt.ioff()
            fig = plt.figure(1)         
            imshow_grid(mg, 'topographic__elevation', grid_units=('m', 'm'), var_name="Elevation (m)", cmap='terrain', allow_colorbar=True)
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/DEM_Image/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=dpi)
            plt.close(fig)
        
        # Slope_Area
        if Slope_Area == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Slope_Area') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Slope_Area')
            plt.ioff()
            fig = plt.figure(2)
            for i, outlet_id in enumerate(prf.data_structure):
                for j, segment_id in enumerate(prf.data_structure[outlet_id]):
                    if j == 0:
                        label = "channel {i}".format(i=i + 1)
                    else:
                        label = '_nolegend_'
                    segment = prf.data_structure[outlet_id][segment_id]
                    profile_ids = segment["ids"]
                    color = segment["color"]
                    plt.loglog(mg.at_node["drainage_area"][profile_ids][2 : -1], mg.at_node["topographic__steepest_slope"][profile_ids][2 : -1], '.', color=color, label=label)
            plt.legend(loc="lower left")
            plt.xlabel("Drainage Area (m^2)")
            plt.ylabel("Channel Slope [m/m]")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.grid(linestyle='--')
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Slope_Area/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=dpi)
            plt.close(fig)
        
        # Channel_Profile
        if Channel_Profile == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Channel_Profile') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Channel_Profile')
            plt.ioff()
            fig = plt.figure(3)
            prf.plot_profiles(field='topographic__elevation', xlabel='Distance Upstream (m)', ylabel='Elevation (m)')
            title_text = '$Year$='+str(total_time)
            plt.title(title_text)
            plt.grid(linestyle='--')
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Channel_Profile/'+str(total_time)+'.'+Export_format,  format=Export_format,  dpi=dpi)
            plt.close(fig)
        
        # Channel_Map
        if Channel_Map == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Channel_Map') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Channel_Map')
            plt.ioff()
            fig = plt.figure(4)
            prf.plot_profiles_in_map_view()
            title_text = '$Year$='+str(total_time)
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Channel_Map/'+str(total_time)+'.'+Export_format,  format=Export_format,  dpi=dpi)
            plt.close(fig)
        
        # Ksn_Profile
        if Ksn_Profile == True:
            prf.run_one_step()
            #sf.calculate_steepnesses()
            if path.exists(str(Directory)+'/TerrainSandbox/Ksn_Profile') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Ksn_Profile')
            plt.ioff()
            fig = plt.figure(5)
            for i, outlet_id in enumerate(prf.data_structure):
                for j, segment_id in enumerate(prf.data_structure[outlet_id]):
                    if j == 0:
                        label = "channel {i}".format(i=i + 1)
                    else:
                        label = '_nolegend_'
                    segment = prf.data_structure[outlet_id][segment_id]
                    profile_ids = segment["ids"]
                    distance_upstream = segment["distances"]
                    color = segment["color"]
                    plt.plot(distance_upstream[2 : -1], mg.at_node["channel__steepness_index"][profile_ids][2 : -1], 'x', color=color, label=label)
            plt.xlabel("Distance Upstream (m)")
            plt.ylabel("Steepness Index")
            plt.legend(loc="upper left")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.grid(linestyle='--')
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Ksn_Profile/'+str(total_time)+'.'+Export_format,  format=Export_format,  dpi=dpi)
            plt.close(fig)
        
        # Ksn_Map
        if Ksn_Map == True:
            prf.run_one_step()
            sf.calculate_steepnesses()
            if path.exists(str(Directory)+'/TerrainSandbox/Ksn_Map') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Ksn_Map')
            plt.ioff()
            fig = plt.figure(6)
            #imshow_grid(mg, "channel__steepness_index", grid_units=("m", "m"), var_name="Steepness Index", cmap="jet")
            Ksn = np.ones(np.size(mg.nodes)) * mg.at_node["channel__steepness_index"]
            Ksn[np.where(mg.node_x >= ((dxy * number_of_columns) - (dxy * 2)))] = 0
            Ksn[np.where(mg.node_y >= ((dxy * number_of_rows) - (dxy * 2)))] = 0
            Ksn[np.where(mg.node_x <= dxy * 2)] = 0
            Ksn[np.where(mg.node_y <= dxy * 2)] = 0
            imshow_grid(mg, Ksn, grid_units=("m", "m"), var_name="Steepness Index", cmap="jet")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Ksn_Map/'+str(total_time)+'.'+Export_format,  format=Export_format,  dpi=dpi)
            plt.close(fig)
        
        # Chi_Profile
        if Chi_Profile == True:
            prf.run_one_step()
            cf.calculate_chi()
            if path.exists(str(Directory)+'/TerrainSandbox/Chi_Profile') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Chi_Profile')
            plt.ioff()
            fig = plt.figure(7)
            for i, outlet_id in enumerate(prf.data_structure):
                for j, segment_id in enumerate(prf.data_structure[outlet_id]):
                    if j == 0:
                        label = "channel {i}".format(i=i + 1)
                    else:
                        label = '_nolegend_'
                    segment = prf.data_structure[outlet_id][segment_id]
                    profile_ids = segment["ids"]
                    color = segment["color"]
                    plt.plot(mg.at_node["channel__chi_index"][profile_ids], mg.at_node["topographic__elevation"][profile_ids], color=color, label=label)
            plt.xlabel("Chi Index (m)")
            plt.ylabel("Elevation (m)")
            plt.legend(loc="lower right")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.grid(linestyle='--')
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Chi_Profile/'+str(total_time)+'.'+Export_format,  format=Export_format,  dpi=dpi)
            plt.close(fig)
        
        # Chi_Map
        if Chi_Map == True:
            prf.run_one_step()
            cf.calculate_chi()
            if path.exists(str(Directory)+'/TerrainSandbox/Chi_Map') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Chi_Map')
            plt.ioff()
            fig = plt.figure(8)
            imshow_grid(mg, "channel__chi_index", grid_units=("m", "m"), var_name="Chi Index", cmap="jet")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Chi_Map/'+str(total_time)+'.'+Export_format,  format=Export_format,  dpi=dpi)
            plt.close(fig)
        
        # Erosion_Rate_Profile
        if Erosion_Rate_Profile == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Erosion_Rate_Profile') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Erosion_Rate_Profile')
            plt.ioff()
            fig = plt.figure(3)
            prf.plot_profiles(field='erosion_rate', xlabel='Distance Upstream (m)', ylabel='Erosion Rate (m/yr)')
            title_text = '$Year$='+str(total_time)
            plt.title(title_text)
            plt.grid(linestyle='--')
            axes = plt.gca()
            axes.set_xlim([dxy,None])
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Erosion_Rate_Profile/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=dpi)
            plt.close(fig)
        
        # Erosion_Rate_Map
        if Erosion_Rate_Map == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Erosion_Rate_Map') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Erosion_Rate_Map')
            plt.ioff()
            fig = plt.figure(8)
            imshow_grid(mg, "erosion_rate", grid_units=("m", "m"), var_name="Erosion Rate (m/yr)", cmap="magma")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Erosion_Rate_Map/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=dpi)
            plt.close(fig)
        
        # DZ_DT_Profile
        if DZ_DT_Profile == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/DZDT_Profile') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/DZDT_Profile')
            plt.ioff()
            fig = plt.figure(3)
            prf.plot_profiles(field='elevational_change', xlabel='Distance Upstream (m)', ylabel='Rate of Elevational Change (m/yr)')
            title_text = '$Year$='+str(total_time)
            plt.title(title_text)
            plt.grid(linestyle='--')
            axes = plt.gca()
            axes.set_xlim([dxy,None])
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/DZDT_Profile/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=dpi)
            plt.close(fig)
        
        # DZ_DT_Map
        if DZ_DT_Map == True:
            #prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/DZDT_Map') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/DZDT_Map')
            plt.ioff()
            fig = plt.figure(8)
            imshow_grid(mg, "elevational_change", grid_units=("m", "m"), var_name="Rate of Elevational Change (m/yr)", cmap="seismic_r", symmetric_cbar = True)
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/DZDT_Map/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=dpi)
            plt.close(fig)
        
        # Ksp_Profile
        if Ksp_Profile == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Ksp_Profile') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Ksp_Profile')
            plt.ioff()
            fig = plt.figure(117)
            if Spatially_Uniform == True or Spatially_Zoned == True:
                prf.plot_profiles(field=Ksp, xlabel='Distance Upstream (m)', ylabel='Ksp')
            if Tilted_Rocks == True or Folded_Rocks == True:
                prf.plot_profiles(field='K_sp', xlabel='Distance Upstream (m)', ylabel='Ksp')
            title_text = '$Year$='+str(total_time) 
            plt.title(title_text)
            plt.grid(linestyle='--')
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Ksp_Profile/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=dpi)
            plt.close(fig)
        
        # Ksp_Map
        if Ksp_Map == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Ksp_Map') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Ksp_Map')
            plt.ioff()
            fig = plt.figure(8)
            if Spatially_Uniform == True or Spatially_Zoned == True:
                imshow_grid(mg, Ksp, grid_units=("m", "m"), var_name="Ksp", cmap="jet")
            if Tilted_Rocks == True or Folded_Rocks == True:
                imshow_grid(mg, 'K_sp', grid_units=("m", "m"), var_name="Ksp", cmap="jet")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Ksp_Map/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=dpi)
            plt.close(fig)  
        
        # Terrain_3D
        if Terrain_3D == True:
            if path.exists(str(Directory)+'/TerrainSandbox/Terrain_3D') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Terrain_3D')
            plt.ioff()
            fig = plt.figure(9)
            ax = plt.axes(projection ='3d') 
            Z = mg.at_node['topographic__elevation'].reshape(mg.shape)
            if Color_Scaling_Updated == True:
                color = cm.terrain((Z-Z.min())/(Z.max()-Z.min()))                  # USE THIS LINE IF YOU WANT THE COLOR SCALING TO UPDATE EACH TIMESTEP ACCORDING TO THE MAX ELEVATION 
            if Color_Scaling_Prescribed == True:
                color = cm.terrain((Z-Z.min())/(Max_color_scale - Z.min()))                     # USE THIS LINE IF YOU WANT THE COLOR SCALING TO STAY CONSTANT THROUGH THE RUN. ENTER DESIRED MAX COLOR ELEVATION WHERE Z.max() SHOULD BE IN THE LINE ABOVE
            if Color_Scaling_Automatic == True:
                color = cm.terrain((Z-Z.min())/(ceiling-Z.min()))                   # USE THIS LINE IF YOU WANT TO AUTOSCALE THE RUN USING "ceiling". NOTE THAT IF THE Ksp OR U PARAMETERS ARE CHANGED MID-RUN, THE SCALING WILL ALSO CHANGE.
            surf = ax.plot_surface(mg.node_x.reshape(mg.shape), mg.node_y.reshape(mg.shape), Z, rstride=1, cstride=1, facecolors=color, edgecolor='black', linewidth=1, antialiased=True)
            ax.view_init(elev=30, azim=-120)                                    # USE THIS LINE TO SET THE FIGURE EYELEVEL AZIMUTH
            ax.set_xlabel('meters')
            ax.set_ylabel('meters')
            ax.set_zlabel('Elevation')
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            if Elevation_Exaggeration_Prescribed == True:
                ax.set_zlim(0, dxy * number_of_rows / Exaggeration_factor)
            if Elevation_Exaggeration_Automatic == True:
                ax.set_zlim(0, ceiling * 10)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Terrain_3D/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=dpi)
            plt.close(fig)
            
        # Mean_Basin_Erosion_Rate
        if Mean_Basin_Erosion_Rate == True:
            if path.exists(str(Directory)+'/TerrainSandbox/Mean_Basin_Erosion_Rate') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Mean_Basin_Erosion_Rate')
            plt.ioff()
            fig = plt.figure()
            x_axis = np.arange(dt, total_time + dt, dt)
            plt.plot(x_axis, mean_basin_erosion_rate, '.--')
            title_text = '$Year$='+str(total_time) 
            plt.title(title_text)
            plt.grid(linestyle='--')
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Mean_Basin_Erosion_Rate/Mean_Basin_Erosion_Rate.' + Export_format,  format = Export_format, dpi=dpi)
            plt.close(fig)   
            
        # Zr_Diff_Plot
        if Zr_Diff == True:
            if path.exists(str(Directory)+'/TerrainSandbox/Zr_Diff_Plot') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Zr_Diff_Plot')
            plt.ioff()
            fig = plt.figure()
            x_axis = np.arange(dt, total_time + dt, dt)
            plt.errorbar(x_axis, Diff['total_zr_diff_mean'], Diff['total_zr_diff_std'], fmt = '.', capsize = 5)
            plt.errorbar(x_axis, Diff['A_zr_diff_mean'], Diff['A_zr_diff_std'], fmt = '.', capsize = 5)
            plt.errorbar(x_axis, Diff['B_zr_diff_mean'], Diff['B_zr_diff_std'], fmt = '.', capsize = 5)
            title_text = '$Year$='+str(total_time) 
            plt.title(title_text)
            plt.grid(linestyle='--')
            plt.tight_layout()
            plt.legend(['Total', 'Zone A', 'Zone B'])
            fig.savefig(str(Directory)+'/TerrainSandbox/Zr_Diff_Plot/Zr_Diff_Plot.' + Export_format,  format = Export_format, dpi=dpi)
            plt.close(fig)          
            
        # Reset Plot_Ticker    
        Plot_Ticker = 0
        
    # Export_DEM    
    if Export_DEM == True:
        if total_time == Export_DEM_Ticker:
            if path.exists(str(Directory)+'/TerrainSandbox') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox')
        if Export_DEM_Ticker == Export_DEM_Interval:
            if path.exists(str(Directory)+'/TerrainSandbox/Export_DEM') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Export_DEM')
            write_esri_ascii(str(Directory)+'/TerrainSandbox/Export_DEM/'+str(total_time)+'.asc', mg, names='topographic__elevation')
            Export_DEM_Ticker = 0
            
print('')
print('Complete!')





#%%

# Print space
print(' ')

# Calculating YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY

# # Set up arrays for analysis
# interval = np.array([100, 250, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000, 100000]) * (dxy**2)
# heads_model = interval * nan
# heads_real = interval * nan
# heads_ratio = interval * nan
# iteration = 0

# for j in interval:

#     # Print status
#     print('j = ', j)
    
#     # Calculate number of channel heads in model   
#     frr = FlowAccumulator(mg, flow_director='D8')
#     frr.run_one_step()
#     prf = ChannelProfiler(mg, number_of_watersheds=1, main_channel_only=False, minimum_channel_threshold=j)
#     prf.run_one_step()         
#     ds = prf.data_structure
#     ds = ds[pour_point_node]
#     heads = []
#     for i in ds.keys():
#         heads = np.append(heads, i[0])
#         heads = np.append(heads, i[1])
#     (unique, counts) = np.unique(heads, return_counts = True)
#     heads_model[iteration] = np.size(np.where(counts == 1)) - 1
        
#     # Calculate number of channel heads in real terrain  
#     frr0 = FlowAccumulator(mg0, flow_director='D8')
#     frr0.run_one_step()   
#     prf0 = ChannelProfiler(mg0, number_of_watersheds=1, main_channel_only=False, minimum_channel_threshold=j)
#     prf0.run_one_step()       
#     ds0 = prf0.data_structure
#     ds0 = ds0[pour_point_node]
#     heads0 = []
#     for i in ds0.keys():
#         heads0 = np.append(heads0, i[0])
#         heads0 = np.append(heads0, i[1])
#     (unique0, counts0) = np.unique(heads0, return_counts = True)
#     heads_real[iteration] = np.size(np.where(counts0 == 1)) - 1   
    
#     # Calculate ratio
#     heads_ratio[iteration] = heads_model[iteration] / heads_real[iteration]
    
#     # Advance iteration
#     iteration += 1

# # Calculate sum
# heads_ratio_sum = np.sum(np.abs(heads_ratio - 1))
# heads_ratio_sum = np.append(heads_ratio_sum, nan)

# Calculating ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
# Plotting YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY

# # Heads_Ratio
# plt.ioff()
# fig = plt.figure()         
# plt.semilogx(interval, heads_ratio, '.-', zorder = 2)
# plt.text(1E7, np.min(heads_ratio), heads_ratio_sum[0])
# plt.semilogx(interval, np.ones(np.size(interval)), 'k--', zorder = 1)
# plt.tight_layout()
# fig.savefig(str(Directory)+'/TerrainSandbox/Analysis/Heads_Ratio.' + Export_format,  format=Export_format, dpi=dpi)
# plt.close(fig)  

# # Zr_Diff
# plt.ioff()
# fig = plt.figure()         
# imshow_grid(mg, zr_diff, grid_units=('m', 'm'), var_name="Elevational Difference (m)", cmap="seismic_r", symmetric_cbar = True)
# plt.tight_layout()
# fig.savefig(str(Directory)+'/TerrainSandbox/Analysis/Zr_Diff. '+ Export_format,  format=Export_format, dpi=dpi)
# plt.close(fig)  

# # DEM_Comparison
# plt.ioff()
# fig = plt.figure()   
# ax1 = plt.subplot(2, 1, 1)      
# imshow_grid(mg0, zr0, grid_units=('m', 'm'), var_name="Elevation (m)", cmap="terrain")
# plt.title('Real')
# ax2 = plt.subplot(2, 1, 2)      
# imshow_grid(mg, zr, grid_units=('m', 'm'), var_name="Elevation (m)", cmap="terrain")      # vmin = np.min(mg0.at_node['topographic__elevation'][mg0.core_nodes]), vmax = np.max(mg0.at_node['topographic__elevation'][mg0.core_nodes])
# plt.title('Model')
# plt.tight_layout()
# fig.savefig(str(Directory)+'/TerrainSandbox/Analysis/DEM_Comparison.' + Export_format,  format=Export_format, dpi=300)
# plt.close(fig)  

# # SA_Comparison
# plt.ioff()
# fig = plt.figure()  
# # 
# ax1 = plt.subplot(2, 1, 1)     
# # 
# frr0 = FlowAccumulator(mg0, flow_director='D8')
# frr0.run_one_step()       
# prf0 = ChannelProfiler(mg0, number_of_watersheds=1, main_channel_only=False, minimum_channel_threshold=(dxy**2) * 2) 
# prf0.run_one_step()
# #
# for i, outlet_id in enumerate(prf0.data_structure):
#     for j, segment_id in enumerate(prf0.data_structure[outlet_id]):
#         if j == 0:
#             label = "channel {i}".format(i=i + 1)
#         else:
#             label = '_nolegend_'
#         segment = prf0.data_structure[outlet_id][segment_id]
#         profile_ids = segment["ids"]
#         color = segment["color"]
#         plt.loglog(mg0.at_node["drainage_area"][profile_ids][2 : -1], mg0.at_node["topographic__steepest_slope"][profile_ids][2 : -1], '.', color=color, label=label)
# plt.legend(loc="lower left")
# plt.xlabel("Drainage Area (m^2)")
# plt.ylabel("Channel Slope [m/m]")
# plt.grid(linestyle='--')
# plt.title('Real')     
# plt.tight_layout()
# #
# ax2 = plt.subplot(2, 1, 2) 
# #
# prf = ChannelProfiler(mg, number_of_watersheds=1, main_channel_only=False, minimum_channel_threshold=(dxy**2) * 2) 
# prf.run_one_step()    
# #
# for i, outlet_id in enumerate(prf.data_structure):
#     for j, segment_id in enumerate(prf.data_structure[outlet_id]):
#         if j == 0:
#             label = "channel {i}".format(i=i + 1)
#         else:
#             label = '_nolegend_'
#         segment = prf.data_structure[outlet_id][segment_id]
#         profile_ids = segment["ids"]
#         color = segment["color"]
#         plt.loglog(mg.at_node["drainage_area"][profile_ids][2 : -1], mg.at_node["topographic__steepest_slope"][profile_ids][2 : -1], '.', color=color, label=label)
# plt.legend(loc="lower left")
# plt.xlabel("Drainage Area (m^2)")
# plt.ylabel("Channel Slope [m/m]")
# plt.grid(linestyle='--')   
# plt.title('Model')
# plt.tight_layout()
# #
# fig.savefig(str(Directory)+'/TerrainSandbox/Analysis/SA_Comparison.' + Export_format,  format=Export_format, dpi=300)
# plt.close(fig) 

# # Channel_Comparison
# plt.ioff()
# fig = plt.figure()   
# ax1 = plt.subplot(2, 1, 1) 
# frr0 = FlowAccumulator(mg0, flow_director='D8')
# frr0.run_one_step()       
# prf0 = ChannelProfiler(mg0, number_of_watersheds=1, main_channel_only=False, minimum_channel_threshold=10000) 
# prf0.run_one_step()
# prf0.plot_profiles_in_map_view()
# plt.title('Real')
# ax2 = plt.subplot(2, 1, 2)  
# prf = ChannelProfiler(mg, number_of_watersheds=1, main_channel_only=False, minimum_channel_threshold=10000) 
# prf.run_one_step()    
# prf.plot_profiles_in_map_view()
# plt.title('Model')
# plt.tight_layout()
# fig.savefig(str(Directory)+'/TerrainSandbox/Analysis/Channel_Comparison.' + Export_format,  format=Export_format, dpi=dpi)
# plt.close(fig)  
                
# Plotting ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^          

# Export CSVs YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY

# if path.exists(str(Directory)+'/TerrainSandbox/CSV') == False:
#     os.mkdir(str(Directory)+'/TerrainSandbox/CSV')              
# for i in Diff:
#     np.savetxt(str(Directory)+'/TerrainSandbox/CSV/' + str(i) + '.csv', Diff[str(i)], delimiter = ",")   
# np.savetxt(str(Directory)+'/TerrainSandbox/CSV/' + 'mean_basin_erosion_rate.csv', mean_basin_erosion_rate, delimiter = ",")   
# # np.savetxt(str(Directory)+'/TerrainSandbox/CSV/' + 'interval.csv', interval, delimiter = ",")   
# # np.savetxt(str(Directory)+'/TerrainSandbox/CSV/' + 'heads_model.csv', heads_model, delimiter = ",")   
# # np.savetxt(str(Directory)+'/TerrainSandbox/CSV/' + 'heads_real.csv', heads_real, delimiter = ",")   
# # np.savetxt(str(Directory)+'/TerrainSandbox/CSV/' + 'heads_ratio.csv', heads_ratio, delimiter = ",")   
# # np.savetxt(str(Directory)+'/TerrainSandbox/CSV/' + 'heads_ratio_sum.csv', heads_ratio_sum, delimiter = ",")   
 
       
# # Export CSVs ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^            

# if path.exists(str(Directory)+'/TerrainSandbox/Export_DEM') == False:
#     os.mkdir(str(Directory)+'/TerrainSandbox/Export_DEM')             
# write_esri_ascii(str(Directory)+'/TerrainSandbox/Export_DEM/'+str(total_time)+'.asc', mg, names='topographic__elevation')
            
              
              
              