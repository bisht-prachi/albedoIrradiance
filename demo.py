# -*- coding: utf-8 -*-
"""
06-05-2024. Prachi Bisht.
demo.
albedo_irradiance.py.

requirements: numpy, pandas, itertools, astropy, plotly(optional)


function: 
    getEarthAlbedodf(filename):
    input:
        1). filename: albedo refelctivity file from 
                      https://neo.gsfc.nasa.gov/view.php?datasetId=MCD43C3_E_BSA&date=2023-12-01
                      can upgrade to higher resolution
    output:
        1). global dataframe containing:
           ['lat', 'lon', 'cell_area', 'albedo']
           for earth on the day 'filename' was collected for

function: 
    getIrradianceAtSat(at_time, sc_x_pos, sc_y_pos, sc_z_pos, filename):
    input:
        1). observation_time: time eg. datetime object e.g. "23-12-2023  00:00:13"
        2). sc_x_pos: spaccraft x poisition (ECEF)
        3). sc_y_pos: spaccraft y poisition (ECEF)
        4). sc_z_pos: spaccraft z poisition (ECEF)
    output:
        1). dataframe containing:
               ['lat', 'lon', 'cell_area', 'albedo', 'earth_radius_at_lat',
              'dot_prod_with_sat', 'sunlit_flag', 'dot_prod_with_sun',
              'dot_prod_with_panel', 'irradiance']
        all within the sc field of view
        
        2). albedo irradiance on a sun-pointing LEO S/C in W/m^2
"""

# import module
import albedoIrradiance as arad
import pandas as pd

# test satellite locations and times
filename = "MCD43C3_M_BSA_2023-12-01_rgb_360x180.SS"
at_time = "2023-12-23  00:00:13"
sc_x_pos, sc_y_pos, sc_z_pos = (237.7391929, 6557.207059, 2746.6659)

# filename = "MCD43C3_M_BSA_2023-03-01_rgb_360x180.SS"
# at_time = "2023-03-20 00:51:00"
# sc_x_pos, sc_y_pos, sc_z_pos = (1207.007071, 961.485658,6935.483990	)

# filename = "MCD43C3_M_BSA_2023-03-01_rgb_360x180.SS"
# at_time = "2023-03-23 01:46:30"
# sc_x_pos, sc_y_pos, sc_z_pos = (-1383.709243, 1007.263423, 6896.115112)

# Initialize func getEarthAlbedodf() with filename to get earth grid with albedo
arad.getEarthAlbedodf(filename)

observation_time = pd.to_datetime(at_time, format="%Y-%m-%d  %H:%M:%S")

# Get irradiance
geo_dataframe, irradiance = arad.getIrradianceAtSat(
    observation_time, sc_x_pos, sc_y_pos, sc_z_pos
)

location = f"({sc_x_pos}, {sc_y_pos}, {sc_z_pos}) (km)"
print(
    f"This is the albedo irradiance at satellite location {location} at time {at_time}:\n{irradiance} W/m^2"
)

# use following if plotly installed
arad.getFOVGeoPlot(geo_dataframe,location, at_time)
arad.getSunlitGeoPlot(geo_dataframe, location, at_time)
arad.getIrradianceGeoPlot(geo_dataframe, location, at_time)
