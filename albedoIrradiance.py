#!/usr/bin/env python
# coding: utf-8

"""
06-05-2024. Prachi Bisht.
albedoIrradiance.py.
"""

# Importing necessary libraries
import numpy as np
import pandas as pd
from itertools import product
import os
from astropy.time import Time
from astropy.coordinates import (
    get_sun,
    AltAz,
    CartesianRepresentation,
    SphericalRepresentation,
    EarthLocation,
)
import astropy.units as u
from astropy import constants as const
import plotly.express as px  # optional


# Getting the directory containing albedo dataset wrt to current directory

## Global constants
am0_intensity = 1367  # W / m^2
earth_mean_radius = const.R_earth.to(u.km)
earth_equatorial_radius = 6378.137 * u.km
earth_polar_radius = 6356.7523142 * u.km


## Functions


def getEarthAlbedodf(albFilename, inFilename):
    # df_Earth is initialized globally. useful if getIrradianceAtSat() is used in loop. need not
    # initialize df_Earth repeatedly
    global df_earth
    
    file_path = os.path.join(os.getcwd(), albFilename)
    df_read = pd.read_csv(file_path, index_col=None, skiprows=8, delim_whitespace=True)


    lon_array = np.linspace(0.5, 359.5, 360)  
    # lon_array = (lon_array + 180) % 360 - 180
    
    lat_array = np.linspace(-89.5, 89.5, 180)
    
    lat_lon_pairs = list(product(lat_array, lon_array))
    
    df_earth = pd.DataFrame(lat_lon_pairs, columns=["lat", "lon"])
    df_earth["albedo"] = df_read['Data']
    
    df_earth["albedo"] = df_earth["albedo"].where(df_earth["albedo"] >= 0)
    
    lat_resolution = abs(lat_array[0] - lat_array[1])
    lon_resolution = abs(lon_array[0] - lon_array[1])
    prefactor = np.radians(lon_resolution)
    
    differential_area = prefactor * np.abs(
        np.cos(np.radians(lat_array - 0.5 * lat_resolution)) -
        np.cos(np.radians(lat_array + 0.5 * lat_resolution))
    )
    
    earth_shape = (len(lat_array), len(lon_array))
    
    cell_area_matrix = np.tile(differential_area[:, np.newaxis], (1, earth_shape[1]))
    df_earth["cell_area"] = cell_area_matrix.flatten() * earth_mean_radius**2

    file_path = os.path.join(os.getcwd(), inFilename)
    df_earth["inSolarFlux"] = pd.read_csv(file_path, index_col=None, skiprows=8, delim_whitespace=True, usecols=["Data"])

    return df_earth


def getRadiusAtLat(df_earth):
    # Calculate Earth radius at each latitude
    f1 = (
        (earth_equatorial_radius**2) * np.cos(np.radians(df_earth["lat"].values))
    ) ** 2
    f2 = ((earth_polar_radius**2) * np.sin(np.radians(df_earth["lat"].values))) ** 2
    f3 = (earth_equatorial_radius * np.cos(np.radians(df_earth["lat"].values))) ** 2
    f4 = (earth_polar_radius * np.sin(np.radians(df_earth["lat"].values))) ** 2

    radius_at_lat = np.sqrt((f1 + f2) / (f3 + f4))
    df_earth["earth_radius_at_lat"] = radius_at_lat

    return df_earth


def getSatVectorECEF(sc_lat, sc_lon, sc_alt):
    # Calculate satellite vector in Cartesian representation
    sat_vector_ecef = SphericalRepresentation(
        lat=sc_lat * u.deg,
        lon=sc_lon * u.deg,
        distance=earth_mean_radius + sc_alt * u.km,
    ).to_cartesian()
    return sat_vector_ecef


def getSatFovdf(df_earth, sat_vector_ecef):
    # Calculate Field of View (FOV) of the satellite
    fov_limit = earth_mean_radius / sat_vector_ecef.norm()

    cell_vectors_cartesian = SphericalRepresentation(
        lat=df_earth["lat"].values * u.deg,
        lon=df_earth["lon"].values * u.deg,
        distance=earth_mean_radius * u.km,
    ).to_cartesian()
    cell_norm_vectors = cell_vectors_cartesian / cell_vectors_cartesian.norm()

    sat_norm_ecef = sat_vector_ecef / sat_vector_ecef.norm()

    df_earth["dot_prod_with_sat"] = sat_norm_ecef.dot(cell_norm_vectors)
    df_earth["satfov_flag"] = df_earth["dot_prod_with_sat"] >= fov_limit
    df = df_earth[df_earth["dot_prod_with_sat"] >= fov_limit].copy()

    return df


def getFOVGeoPlot(location, at_time):
    fig = px.scatter_geo(
        df_earth,
        lat="lat",
        lon="lon",
        color="satfov_flag",
        opacity=0.3,
        height=600,
        color_continuous_scale="plasma_r",
        labels={"satfov": ""},
    )
    fig.update_traces(marker=dict(size=5))
    fig.update_geos(resolution=110)
    fig.update_layout(
        title_text="<b>satellite field-of-view </b><br>",
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
        ),
    )
    formatted_time = at_time.replace(":", "-")
    # geo visulaization of the field-of-view of satellite saved as html file
    fig.write_html(f"satfov_{formatted_time}.html")
    fig.write_image(f"satfov_{formatted_time}.png")
    fig.show()


def getFOVElementVectorsECEF(df):
    # Calculate element vectors in ECEF coordinates
    element_vectors_ecef = SphericalRepresentation(
        lat=df["lat"].values * u.deg,
        lon=df["lon"].values * u.deg,
        distance=earth_mean_radius,
    ).to_cartesian()

    return element_vectors_ecef


def getFOVElementVectorsECEF(df):
    # Calculate element vectors in ECEF coordinates
    element_vectors_ecef = SphericalRepresentation(
        lat=df["lat"].values * u.deg,
        lon=df["lon"].values * u.deg,
        distance=earth_mean_radius,
    ).to_cartesian()

    return element_vectors_ecef


def getFOVDotProductwithSun(df, observation_time):
    # Get sunlit_flag and calculate dot product with the sun
    locations = EarthLocation(
        lat=df["lat"].values * u.deg, lon=df["lon"].values * u.deg, height=0.0
    )

    sun_gcrs = get_sun(observation_time)
    sun_altaz = sun_gcrs.transform_to(AltAz(location=locations))
    df["sunlit_flag"] = sun_altaz.alt.deg > 0

    # sun_altaz_cartesian = sun_altaz.represent_as(CartesianRepresentation)
    # sun_norm_altaz = sun_altaz_cartesian / sun_altaz_cartesian.norm()

    sun_vector_ecef = sun_gcrs.transform_to("itrs").cartesian
    sun_norm_ecef = sun_vector_ecef / sun_vector_ecef.norm()

    element_vectors_ecef = getFOVElementVectorsECEF(df)
    element_norm_ecef = element_vectors_ecef / element_vectors_ecef.norm()
    
    # df["dot_prod_with_sun"] = sun_norm_altaz.dot(element_norm_ecef)
    df["dot_prod_with_sun"] = sun_norm_ecef.dot(element_norm_ecef)

    sun_vector_ecef = sun_gcrs.transform_to("itrs").cartesian

    return sun_vector_ecef, df


def getFOVDotProductwithPanel(df, sat_vector_ecef, sun_vector_ecef):
    # Calculate dot product with the panel
    sat_norm_ecef = sat_vector_ecef / sat_vector_ecef.norm()
    sun_norm_ecef = sun_vector_ecef / sun_vector_ecef.norm()
    
    panel_vector_ecef = (
         sat_norm_ecef - sun_norm_ecef
    )
    # panel_vector_ecef = sat_vector_ecef      #if sc is earth-pointing
    panel_norm_ecef = panel_vector_ecef / panel_vector_ecef.norm()

    element_vectors_ecef = getFOVElementVectorsECEF(df)
    element_norm_ecef = element_vectors_ecef / element_vectors_ecef.norm()

    df["dot_prod_sat_with_panel"] = panel_norm_ecef.dot(sat_norm_ecef)
    df["dot_prod_with_panel"] = panel_norm_ecef.dot(element_norm_ecef)

    return panel_vector_ecef, df


def getDistToSat(df, sat_vector_ecef):
    # Calculate the distance between the elements and satellite
    element_vectors_ecef = getFOVElementVectorsECEF(df)
    df["dist_to_sat"] = (sat_vector_ecef - element_vectors_ecef).norm()

    return df


def getIrradianceAtSat(at_time, sc_lat, sc_lon, sc_alt):
    # Main function to calculate irradiance at the satellite
    observation_time = Time(at_time)  # , '%d-%m-%Y  %H:%M:%S')
    sat_vector_ecef = getSatVectorECEF(sc_lat, sc_lon, sc_alt)

    df = getSatFovdf(df_earth, sat_vector_ecef)

    sun_vector_ecef, df = getFOVDotProductwithSun(df, observation_time)
    panel_vector_ecef, df = getFOVDotProductwithPanel(
        df, sat_vector_ecef, sun_vector_ecef
    )

    df = getDistToSat(df, sat_vector_ecef)
    # Calculate irradiance
    df["irradiance"] = (
        1
        * (
            df["sunlit_flag"]
            * df["albedo"]
            * df["cell_area"]
            * df["dot_prod_with_sun"]
            * df["dot_prod_with_sat"]
            * df_earth["inSolarFlux"]
        )
        / (np.pi * df["dist_to_sat"] ** 2)
    )

    # df["irradiance"] = df["irradiance"].clip(0, None)

    irradiance = df["irradiance"].sum()

    return irradiance, df

def getAlbedoGeoPlot(df, location, at_time):
    fig = px.scatter_geo(
        df,
        lat="lat",
        lon="lon",
        color="albedo",
        height=600,
        color_continuous_scale="Blues_r",
        opacity=0.5,
        labels={"sunlit": ""},
    )
    fig.update_traces(marker=dict(size=5))
    fig.update_geos(resolution=110)
    fig.update_layout(
        title_text="<b>albedo</b><br>",
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
        ),
    )
    formatted_time = at_time.replace(":", "-")
    # geo visulaization of sunlit area the field-of-view of satellite saved as html file
    fig.write_html(f"albedo_{formatted_time}.html")
    fig.write_image(f"albedo_{formatted_time}.png")
    fig.show()
    

def getSunlitGeoPlot(df, location, at_time):
    fig = px.scatter_geo(
        df,
        lat="lat",
        lon="lon",
        color="sunlit_flag",
        height=600,
        color_continuous_scale="Blues_r",
        opacity=0.5,
        labels={"sunlit": ""},
    )
    fig.update_traces(marker=dict(size=5))
    fig.update_geos(resolution=110)
    fig.update_layout(
        title_text="<b>sunlit area</b><br>",
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
        ),
    )
    formatted_time = at_time.replace(":", "-")
    # geo visulaization of sunlit area the field-of-view of satellite saved as html file
    fig.write_html(f"sunlit_{formatted_time}.html")
    fig.write_image(f"sunlit_{formatted_time}.png")
    fig.show()
    


def getIrradianceGeoPlot(df, location, at_time):
    fig = px.scatter_geo(
        df,
        lat="lat",
        lon="lon",
        color="irradiance",
        opacity=0.3,
        height=600,
        color_continuous_scale="Blues_r",
        range_color=[df["irradiance"].min(), df["irradiance"].max()],
        labels={"irradiance": ""},
    )

    fig.update_traces(marker=dict(size=5))
    fig.update_geos(resolution=110)
    fig.update_layout(
        title_text=f"<b>irradiance (W/m^2) from satellite FOV <br> satellite at {location} on {at_time}</b><br>",
        title_x=0.5,
        title_y=0.97,
        font_family="Arial",
        font_size=16,
        geo=dict(bgcolor="cornflowerblue", landcolor="forestgreen"),
        coloraxis_colorbar=dict(
            orientation="v",
            len=0.65,
            xanchor="right",
            x=1.08,
            yanchor="bottom",
            y=0.15,
            thickness=25,
            bgcolor="white",
        ),
    )

    formatted_time = at_time.replace(":", "-")
    # geo visulaization of irradiance from the field-of-view of satellite saved as html file
    fig.write_html(f"irradiance_{formatted_time}.html")
    fig.write_image(f"irradiance_{formatted_time}.png")
    fig.show()



def main():
    # Main function to execute the code

    # test satellite locations and times
    month = "December"
    year = "2023"
    albFilename = f"CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed4.1_Observed_TOA_Albedo-All-sky_{month}-{year}{month}-{year}.txt"#"MCD43C3_M_BSA_2023-12-01_rgb_360x180.SS"
    inFilename = f"CERES_EBAF-TOA_Ed4.2_Incoming_Solar_Flux_{month}-{year}{month}-{year}.txt"
    
    at_time = "2023-12-23  00:00:13"
    sc_lat, sc_lon, sc_alt = (22, 88, 760)

    # Initialize func getEarthAlbedodf() with filename to get earth grid with albedo
    getEarthAlbedodf(albFilename, inFilename)

    observation_time = pd.to_datetime(at_time, format="%Y-%m-%d  %H:%M:%S")

    # Get irradiance
    irradiance, geo_dataframe = getIrradianceAtSat(
        observation_time, sc_lat, sc_lon, sc_alt
    )

    location = f"({sc_lat}\u00b0, {sc_lon}\u00b0, {sc_alt} km) "
    print(
        f"This is the albedo irradiance at satellite location {location} at time {at_time}:\n{irradiance} W/m^2"
    )


if __name__ == "__main__":
    main()
