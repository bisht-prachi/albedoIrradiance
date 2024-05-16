# Albedo Irradiance Calculator

Author: Prachi Bisht
Date: 06-05-2024

## Overview
albedoIrradiance.py script calculates the irradiance due to earthshine/albedo on a LEO spacecraft 
at a given time and (latitude, longitude, altitude)
using reflectivity dataset from Moderate Resolution Imaging Spectroradiometer (MODIS) aboard Terra Satellite. 
Use case: quantify the earthshine irradiance on photovoltaics like solar panels/ sun sensors aboard a LEO spacecraft 

See 'demo.py' for guidance.


## Function

getEarthAlbedodf(filename)

	Input:
	1). filename: albedo refelctivity file from 
	https://neo.gsfc.nasa.gov/view.php?datasetId=MCD43C3_E_BSA&date=2023-12-01
	select 'CSV for Excel' file type with 360x180 resolution.
	can upgrade to higher resolution
 
	Output:
	1). global dataframe containing:
	['lat', 'lon', 'cell_area', 'albedo']
	of earth for the month 'filename' was collected for

## Function

getIrradianceAtSat(at_time, sc_lat, sc_lon, sc_alt)

	Input:
	1). at_time: time eg. "23-12-2023  00:00:13"
	2). sc_lat: spaccraft latitude (degrees)
        3). sc_lon: spaccraft longitude (degrees)
        4). sc_alt: spaccraft altitude (km)

	Output:
	1). dataframe containing:
	       ['lat', 'lon', 'cell_area', 'albedo',
	      'dot_prod_with_sun', 'sunlit_flag', 
	      'dot_prod_with_sat', 'irradiance']
	    all within the sc field of view
	
	2). albedo irradiance on a sun-pointing sc in W/m^2

## Requirements

- numpy==1.21.5
- pandas==1.3.3
- astropy==4.3.1

You can install these dependencies using pip:

```terminal
pip install -r requirements.txt
```

## Usage

1. Place your albedo data CSV file in the subdirectory 'albedo_dataset'.
2. Update the `filename` (albedo data CSV file), `sc_lat`, `sc_lon`, `sc_alt`, and `at_time` variables in the `main` function with your desired values.
3. Run the script:

```terminal
python albedoIrradiance.py
```

4. The script will output the irradiance at the specified location and time in Watts per square meter (W/m^2).

## Example

Location: (22°, 88°, 740km)
Time: 23-12-2023 00:00:13
Irradiance: [irradiance] W/m^2
![satfov_2023-12-23  00-00-13](https://github.com/bisht-prachi/albedoIrradiance/assets/103419553/23e7d73e-225b-4e5d-9cfe-1c7da470899a)
![sunlit_2023-12-23  00-00-13](https://github.com/bisht-prachi/albedoIrradiance/assets/103419553/63ab9fc4-2657-4108-b7a1-173e865f3ac5)
![irradiance_2023-12-23 00-00-13](https://github.com/bisht-prachi/albedoIrradiance/assets/103419553/cb4c1475-b1bd-48a7-be5c-56226c25d4a7)


