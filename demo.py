# -*- coding: utf-8 -*-
"""
06-05-2024. Prachi Bisht.
demo.
albedo_irradiance.py.

requirements: numpy, pandas, itertools, astropy, plotly(optional)


function: 
    getEarthAlbedodf(filename):
    input:
        1). albFilename: albedo file from 
                      https://ceres.larc.nasa.gov/data/#syn1deg-level-3 (SYN1deg)
                      for a specified date

    
        2). inFilename: incoming solar flux file from 
                      https://ceres.larc.nasa.gov/data/#syn1deg-level-3 (SYN1deg)
                      for a specified date

    output:
        1). global dataframe containing: ['lat', 'lon', 'cell_area', 'albedo']
                                    	of earth for the date/month 'inFilename'/'albFilename' 
                                        was collected for

function: 
    getIrradianceAtSat(at_time, sc_x_pos, sc_y_pos, sc_z_pos, filename):
    input:
        1). observation_time: time eg. datetime object e.g. "23-12-2023  00:00:13"
        2). sc_lat: spacecraft latitude (degrees)
        3). sc_lon: spacecraft longitude (degrees)
        4). sc_alt: spacecraft altitude (km)
    output:
        1). dataframe containing:
               ['lat', 'lon', 'cell_area', 'albedo', 'earth_radius_at_lat',
              'dot_prod_with_sat', 'sunlit_flag', 'dot_prod_with_sun',
              'dot_prod_with_panel', 'irradiance']
        all within the S/C field of view
        
        2). albedo irradiance on LEO S/C in W/m^2
"""

# import module
import albedoIrradiance as arad
import pandas as pd

# test satellite locations and times
month = "December"
year = "2023"
albFilename = f"CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed4.1_Observed_TOA_Albedo-All-sky_{month}-{year}{month}-{year}.txt"#"MCD43C3_M_BSA_2023-12-01_rgb_360x180.SS"
inFilename = f"CERES_EBAF-TOA_Ed4.2_Incoming_Solar_Flux_{month}-{year}{month}-{year}.txt"

at_time = "2023-12-23  00:00:13"
sc_lat, sc_lon, sc_alt = (22, 88, 760)

# Initialize func getEarthAlbedodf() with filename to get earth grid with albedo
arad.getEarthAlbedodf(albFilename, inFilename)

observation_time = pd.to_datetime(at_time, format="%Y-%m-%d  %H:%M:%S")

# Get irradiance and dataframe containing irradiance from the FOV
irradiance, geo_dataframe = arad.getIrradianceAtSat(
    observation_time, sc_lat, sc_lon, sc_alt
)

location = f"({sc_lat}\u00b0, {sc_lon}\u00b0, {sc_alt} km) "
print(
    f"This is the albedo irradiance at satellite location {location} at time {at_time}:\n{irradiance} W/m^2"
)

# use following if plotly installed
arad.getFOVGeoPlot(location, at_time)
arad.getSunlitGeoPlot(geo_dataframe, location, at_time)
arad.getAlbedoGeoPlot(geo_dataframe, location, at_time)
arad.getIrradianceGeoPlot(geo_dataframe, location, at_time)
