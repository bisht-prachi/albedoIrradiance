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
import plotly.express as px  # (optional)

# test satellite location and time
filename = "\\MCD43C3_E_BSA_2023-12-19_rgb_360x180.SS"
at_time = "2023-12-23 20:44:00"
sc_x_pos, sc_y_pos, sc_z_pos = (
    -85.757013,
    -6696.066279,
    -2415.234881,
) 

observation_time = pd.to_datetime(at_time, format="%Y-%m-%d  %H:%M:%S")

# Initialize func getEarthAlbedodf() with filename to get earth grid with albedo

arad.getEarthAlbedodf(filename)

# Get irradiance
geo_dataframe, irradiance = arad.getIrradianceAtSat(
    observation_time, sc_x_pos, sc_y_pos, sc_z_pos
)

location = f"({sc_x_pos}, {sc_y_pos}, {sc_z_pos}) (km)"
print(
    f"This is the albedo irradiance at satellite location {location} at time {at_time}:\n{irradiance} W/m^2"
)


# use following if plotly installed
fig = px.scatter_geo(
    geo_dataframe,
    lat="lat",
    lon="lon",
    color="irradiance",
    opacity=0.4,
    height=600,
    color_continuous_scale="Blues_r",
    range_color=[geo_dataframe["irradiance"].min(), geo_dataframe["irradiance"].max()],
    labels={"irradiance": ""},
)

fig.update_traces(marker=dict(size=5))
fig.update_geos(resolution=110)
fig.update_layout(
    title_text="<b>irradiance W/m^2 </b>",
    title_x=0.5,
    title_y=0.93,
    font_family="Arial",
    font_size=22,
    geo=dict(bgcolor="cornflowerblue", landcolor="forestgreen"),
    coloraxis_colorbar=dict(
        orientation="v",
        len=0.65,
        xanchor="right",
        x=1.08,
        yanchor="bottom",
        y=0.15,
        thickness=35,
        bgcolor="white",
    ),
)

# geo visulaization of irradiance from the field-of-view of satellite saved as html file
fig.write_html("irradiance_23122024.html")

fig.show()
