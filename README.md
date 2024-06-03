# Albedo Irradiance Calculator

Author: Prachi Bisht
Date: 06-05-2024

## Overview
albedoIrradiance.py script calculates the irradiance due to earthshine/albedo on a LEO spacecraft 
at a given time and (latitude, longitude, altitude)
using CERES (Clouds and the Earth’s Radiant Energy System) Data Product SYN1Deg. 

Use case: quantify the earthshine irradiance on photovoltaics like solar panels/ sun sensors aboard a LEO spacecraft 

See 'demo.py' for guidance.


## Function

getEarthAlbedodf(filename)

	Input:
	input:
        1). albFilename: albedo refelctivity file from 
                      https://ceres.larc.nasa.gov/data/#syn1deg-level-3 (SYN1deg)
    
        2). inFilename: incoming solar flux file from 
                      https://ceres.larc.nasa.gov/data/#syn1deg-level-3 (SYN1deg)
 
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
![satfov_2023-12-23  00-00-13](https://github.com/bisht-prachi/albedoIrradiance/assets/103419553/8afa6d9b-d983-43b8-abc8-4bf6ac1a10e0)
![albedo_2023-12-23  00-00-13](https://github.com/bisht-prachi/albedoIrradiance/assets/103419553/8aaf5589-a864-47c2-a607-209dd007be79)
![sunlit_2023-12-23  00-00-13](https://github.com/bisht-prachi/albedoIrradiance/assets/103419553/2383ae89-ebde-445b-9152-42387acd6e19)
![irradiance_2023-12-23  00-00-13](https://github.com/bisht-prachi/albedoIrradiance/assets/103419553/6f3c925f-affb-494e-a222-d5988ffdc237)




