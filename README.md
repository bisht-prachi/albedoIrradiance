# Albedo Irradiance Calculator

Author: Prachi Bisht
Date: 06-05-2024

## Overview
albedoIrradiance.py script calculates the albedo irradiance on a LEO spacecraft
with sun-pointing solar panels at a given location and time using albedo
dataset from Moderate Resolution Imaging Spectroradiometer (MODIS) aboard Terra Satellite. 

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

getIrradianceAtSat(at_time, sc_x_pos, sc_y_pos, sc_z_pos)

	Input:
	1). at_time: time eg. "23-12-2023  00:00:13"
	2). sc_x_pos: spaccraft x poisition (ECEF)
	3). sc_y_pos: spaccraft y poisition (ECEF)
	4). sc_z_pos: spaccraft z poisition (ECEF)

	Output:
	1). dataframe containing:
	       ['lat', 'lon', 'cell_area', 'albedo', 'earth_radius_at_lat',
	      'dot_prod_with_sat', 'sunlit_flag', 'dot_prod_with_sun',
	      'dot_prod_with_panel', 'irradiance']
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
2. Update the `filename` (albedo data CSV file), `sc_x_pos`, `sc_y_pos`, `sc_z_pos`, and `at_time` variables in the `main` function with your desired values.
3. Run the script:

```terminal
python albedoIrradiance.py
```

4. The script will output the irradiance at the specified location and time in Watts per square meter (W/m^2).

## Example

Location: (237.7391929, 6557.207059, 2746.6659)
Time: 23-12-2023 00:00:13
Irradiance: [irradiance value] W/m^2
![satfov_23122024](https://github.com/bisht-prachi/albedoIrradiance/assets/103419553/71a18c45-1a12-4c9d-b7e1-1327da979887)

![sunlit_23122024](https://github.com/bisht-prachi/albedoIrradiance/assets/103419553/ef3e5b81-d698-4410-b6aa-8f48fa3ebac9)

![irradiance_2023-12-23 00-00-13](https://github.com/bisht-prachi/albedoIrradiance/assets/103419553/896dcb2e-f322-47f5-9eff-6d145cf016c3)







