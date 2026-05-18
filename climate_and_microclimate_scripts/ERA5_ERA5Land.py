import os
import pandas as pd
import numpy as np
import xarray as xr
import datetime
import calendar
import glob
import re
import cfgrib
import warnings
import json
from scipy.interpolate import griddata
import shutil
import zipfile
import argparse

# Suppress warnings
warnings.filterwarnings('ignore')

class ERA5CombinedProcessor:
    """
    Complete ERA5 and ERA5-Land data processing workflow with proper grid alignment,
    intermediate file handling, and comprehensive diagnostics.
    """
    
    def __init__(self, base_dir="era5_workflow_new", year=2015, test_mode=False, force_overwrite=False):
        """
        Initialize the ERA5 Combined Processor.
        
        Args:
            base_dir (str): Base directory for all files
            year (int): Year to process
            test_mode (bool): If True, process only 5 days of June for faster testing
            force_overwrite (bool): If True, overwrite existing files and restart from beginning
        """
        # Setup processing parameters
        self.year = year
        self.summer_months = [6, 7, 8]  # June, July, August
        self.test_mode = test_mode
        self.force_overwrite = force_overwrite
        
        # Use a common base directory for all years
        self.base_dir = base_dir
        
        # Create timestamp for this run
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Setup directory structure - use year-specific subdirectories for processing
        self.year_dir = os.path.join(self.base_dir, f"year_{year}")
        self.grib_dir = os.path.join(self.year_dir, "grib_files")
        self.csv_dir = os.path.join(self.year_dir, "csv_files") 
        self.intermediate_dir = os.path.join(self.year_dir, "intermediate")
        self.netcdf_dir = os.path.join(self.year_dir, "netcdf_files")
        self.log_dir = os.path.join(self.base_dir, "logs")
        
        # Define the merged output directory - common for all years
        self.merged_dir = os.path.join(self.base_dir, "merged")
        self.merged_corrected_dir = os.path.join(self.base_dir, "merged_corrected")
        
        # Create directories
        for directory in [self.base_dir, self.year_dir, self.grib_dir, self.csv_dir, 
                          self.intermediate_dir, self.netcdf_dir, self.merged_dir, 
                          self.merged_corrected_dir, self.log_dir]:
            os.makedirs(directory, exist_ok=True)
            
        # Create subdirectories for GRIB files
        self.era5_dir = os.path.join(self.grib_dir, "era5")
        self.era5land_dir = os.path.join(self.grib_dir, "era5-land")
        os.makedirs(self.era5_dir, exist_ok=True)
        os.makedirs(self.era5land_dir, exist_ok=True)
        
        # Resolution information
        self.era5_res = 0.25  # ERA5 resolution in degrees (0.25° × 0.25°)
        self.era5land_res = 0.1  # ERA5-Land resolution in degrees (0.1° × 0.1°)
        
        # Area of interest - Colorado Front Range [N, W, S, E]
        north = 40.5  # Aligned to 0.1° grid
        west = -106.0  # Aligned to 0.1° grid
        south = 39.5  # Aligned to 0.1° grid
        east = -105.0  # Aligned to 0.1° grid
        
        self.area = [north, west, south, east]
        
        # Setup logging
        self.log_file = os.path.join(self.log_dir, f"era5_processor_{year}_{timestamp}.log")
        
        # Initialize CDS API client if the module is available
        try:
            import cdsapi
            self.client = cdsapi.Client()
            self.cdsapi_available = True
        except ImportError:
            self.log("cdsapi module not available. Download functionality disabled.")
            self.cdsapi_available = False
        
        # Log initialization
        self.log(f"ERA5 Combined Processor initialized for year {self.year}")
        self.log(f"Base directory: {self.base_dir}")
        self.log(f"Year-specific directory: {self.year_dir}")
        self.log(f"Merged output directory: {self.merged_dir}")
        self.log(f"Area of interest: {self.area}")
        self.log(f"Test mode: {self.test_mode}")
        
        # Variable mapping from ERA5/ERA5-Land to micro_era5 names
        self.variable_mapping = {
            'temperature_2m': 't2m',
            'dewpoint_2m': 'd2m',
            'surface_pressure': 'sp',
            'wind_u_10m': 'u10',
            'wind_v_10m': 'v10',
            'total_precipitation': 'tp',
            'solar_radiation_downward': 'ssrd',
            'direct_solar_radiation': 'fdir',
            'total_cloud_cover': 'tcc',
            'longwave_radiation_downward': 'msdwlwrf',
            'longwave_radiation_net': 'msnlwrf',
            'avg_sdlwrf': 'msdwlwrf',
            'avg_snlwrf': 'msnlwrf'
        }
        
        # ERA5 variable names in API request
        self.era5_variables = {
            'total_cloud_cover': 'tcc',
            'mean_surface_net_long_wave_radiation_flux': 'msnlwrf',
            'mean_surface_downward_long_wave_radiation_flux': 'msdwlwrf',
            'total_sky_direct_solar_radiation_at_surface': 'fdir'
        }
        
        # ERA5-Land variable names in API request
        self.era5land_variables = {
            '2m_temperature': 't2m',
            '2m_dewpoint_temperature': 'd2m',
            'surface_pressure': 'sp',
            '10m_u_component_of_wind': 'u10',
            '10m_v_component_of_wind': 'v10',
            'total_precipitation': 'tp',
            'surface_solar_radiation_downwards': 'ssrd'
        }
        
        # Define the correct order of variables for micro_era5 compatibility
        self.variable_order = [
            'tp', 'fdir', 'ssrd', 'msnlwrf', 'msdwlwrf', 
            't2m', 'd2m', 'sp', 'u10', 'v10', 'tcc', 'lsm'
        ]
        
        # Variables to interpolate from ERA5 radiation
        self.rad_vars = ['tcc', 'fdir', 'msdwlwrf', 'msnlwrf']
        
        # Variables to keep from ERA5-Land
        self.land_vars = ['sp', 'u10', 'v10', 't2m', 'd2m', 'ssrd', 'tp']
        
        # Default values for filling missing data
        self.default_values = {
            'tcc': 0.5,            # 50% cloud cover
            'fdir': 70000.0,       # Direct radiation
            'msdwlwrf': 300.0,     # Downward longwave
            'msnlwrf': 300.0,      # Net longwave
            'sp': 76000.0,         # Surface pressure for mountains
            'u10': 0.0,            # Wind u component
            'v10': 0.0,            # Wind v component
            't2m': 283.0,          # Temperature (10°C)
            'd2m': 283.0,          # Dewpoint (10°C)
            'ssrd': 100000.0,      # Solar radiation
            'tp': 0.0              # Precipitation
        }
    
    def log(self, message):
        """Log a message to both console and log file"""
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_message = f"[{timestamp}] {message}"
        
        print(log_message)
        
        with open(self.log_file, 'a') as f:
            f.write(log_message + "\n")
    
    def detect_current_step(self):
        """Detect what step the processing is currently at"""
        if self.force_overwrite:
            self.log("Force overwrite enabled - starting from step 1")
            return 1
        
        # Check for corrected files (step 7 complete)
        corrected_nc = os.path.join(self.merged_corrected_dir, f"era5_summer_{self.year}_corrected.nc")
        if os.path.exists(corrected_nc):
            self.log("Step 7 (SSRD and TP correction) already complete")
            return 8  # All done
        
        # Check for final NetCDF (step 5 complete)
        final_nc = os.path.join(self.merged_dir, f"era5_summer_{self.year}.nc")
        if os.path.exists(final_nc):
            self.log("Step 5 (final NetCDF) complete - starting from step 7 (SSRD and TP correction)")
            return 7
        
        # Check for merged data (step 4 complete)
        merged_csv = os.path.join(self.intermediate_dir, "merged_interpolated.csv")
        if os.path.exists(merged_csv):
            self.log("Step 4 (merging) complete - starting from step 5 (final NetCDF)")
            return 5
        
        # Check for CSV analysis (step 3 complete)
        analysis_file = os.path.join(self.intermediate_dir, 'csv_analysis.txt')
        if os.path.exists(analysis_file):
            self.log("Step 3 (CSV analysis) complete - starting from step 4 (merging)")
            return 4
        
        # Check for CSV files (step 2 complete)
        csv_inventory = os.path.join(self.intermediate_dir, 'csv_inventory.json')
        if os.path.exists(csv_inventory):
            with open(csv_inventory, 'r') as f:
                csv_data = json.load(f)
                if csv_data.get('era5_rad') and csv_data.get('era5_land'):
                    # Verify files still exist
                    all_exist = all(os.path.exists(f) for f in csv_data['era5_rad'] + csv_data['era5_land'])
                    if all_exist:
                        self.log("Step 2 (GRIB to CSV) complete - starting from step 3 (analysis)")
                        return 3
        
        # Check for GRIB files (step 1 complete)
        grib_inventory = os.path.join(self.intermediate_dir, 'grib_inventory.json')
        if os.path.exists(grib_inventory):
            with open(grib_inventory, 'r') as f:
                grib_data = json.load(f)
                if grib_data.get('era5') and grib_data.get('era5_land'):
                    # Verify files still exist
                    all_files = list(grib_data['era5'].values()) + list(grib_data['era5_land'].values())
                    all_exist = all(os.path.exists(f) for f in all_files)
                    if all_exist:
                        self.log("Step 1 (GRIB files) complete - starting from step 2 (conversion)")
                        return 2
        
        self.log("Starting from step 1 (download/locate GRIB files)")
        return 1
    
    def correct_tp(self, df):
        """
        Correct Total Precipitation from cumulative m to hourly m.
        
        This handles the fact that tp in ERA5-Land is cumulative since forecast start time
        and resets at midnight (00:00 UTC). We convert it to hourly rates keeping units in m.
        
        Args:
            df (pd.DataFrame): DataFrame with 'valid_time' and 'tp' columns
            
        Returns:
            pd.DataFrame: DataFrame with corrected 'tp_corrected' column
        """
        if 'tp' not in df.columns or 'valid_time' not in df.columns:
            self.log("Warning: Cannot correct tp - missing required columns")
            if 'tp' in df.columns:
                df['tp_corrected'] = df['tp']  # Use original if can't correct
            return df
        
        self.log("Applying tp correction (keeping m units)...")
        
        # Ensure valid_time is datetime
        df['valid_time'] = pd.to_datetime(df['valid_time'])
        
        # Sort by time and location
        df_sorted = df.sort_values(['latitude', 'longitude', 'valid_time']).copy()
        
        # Group by location to handle each grid point separately
        corrected_dfs = []
        
        unique_locations = df_sorted[['latitude', 'longitude']].drop_duplicates()
        total_locations = len(unique_locations)
        
        for idx, (_, location) in enumerate(unique_locations.iterrows()):
            if idx % 100 == 0:  # Log progress
                self.log(f"Processing tp correction for location {idx+1}/{total_locations}")
            
            lat, lon = location['latitude'], location['longitude']
            location_df = df_sorted[
                (df_sorted['latitude'] == lat) & (df_sorted['longitude'] == lon)
            ].copy()
            
            if len(location_df) == 0:
                continue
                
            location_df = location_df.sort_values('valid_time')
            
            # Extract hour
            location_df['hour'] = location_df['valid_time'].dt.hour
            
            # Calculate differences between consecutive tp values
            location_df['tp_diff'] = location_df['tp'].diff().fillna(0)
            
            # Apply correction logic (same as SSRD)
            location_df['tp_hourly'] = 0.0
            
            for i in range(len(location_df)):
                hour = location_df.iloc[i]['hour']
                tp_val = location_df.iloc[i]['tp']
                tp_diff = location_df.iloc[i]['tp_diff']
                
                if hour == 1:  # 01:00 UTC (first hour after midnight reset)
                    location_df.iloc[i, location_df.columns.get_loc('tp_hourly')] = tp_val
                elif tp_diff < 0 and hour != 1:  # Handle resets that aren't at 01:00
                    location_df.iloc[i, location_df.columns.get_loc('tp_hourly')] = 0
                else:
                    location_df.iloc[i, location_df.columns.get_loc('tp_hourly')] = tp_diff
            
            # Keep in m - no conversion needed
            location_df['tp_corrected'] = location_df['tp_hourly']
            
            # Remove intermediate columns
            location_df = location_df.drop(columns=['hour', 'tp_diff', 'tp_hourly'])
            
            corrected_dfs.append(location_df)
        
        # Combine all corrected data
        result_df = pd.concat(corrected_dfs, ignore_index=True)
        
        # Ensure we maintain the original order
        result_df = result_df.sort_values(['valid_time', 'latitude', 'longitude'])
        
        self.log(f"tp correction complete. Original range: [{df['tp'].min():.6f}, {df['tp'].max():.6f}] m")
        if 'tp_corrected' in result_df.columns:
            self.log(f"Corrected range: [{result_df['tp_corrected'].min():.6f}, {result_df['tp_corrected'].max():.6f}] m")
        
        return result_df
    
    def correct_ssrd_joules(self, df):
        """
        Correct SSRD (Solar Surface Radiation Downwards) from cumulative J/m² to hourly J/m².
        
        This handles the fact that SSRD in ERA5 is cumulative since forecast start time
        and resets at midnight (00:00 UTC). We convert it to hourly rates but keep units in J/m².
        
        Args:
            df (pd.DataFrame): DataFrame with 'valid_time' and 'ssrd' columns
            
        Returns:
            pd.DataFrame: DataFrame with corrected 'ssrd_corrected' column
        """
        if 'ssrd' not in df.columns or 'valid_time' not in df.columns:
            self.log("Warning: Cannot correct SSRD - missing required columns")
            if 'ssrd' in df.columns:
                df['ssrd_corrected'] = df['ssrd']  # Use original if can't correct
            return df
        
        self.log("Applying SSRD correction (keeping J/m² units)...")
        
        # Ensure valid_time is datetime
        df['valid_time'] = pd.to_datetime(df['valid_time'])
        
        # Sort by time and location
        df_sorted = df.sort_values(['latitude', 'longitude', 'valid_time']).copy()
        
        # Group by location to handle each grid point separately
        corrected_dfs = []
        
        unique_locations = df_sorted[['latitude', 'longitude']].drop_duplicates()
        total_locations = len(unique_locations)
        
        for idx, (_, location) in enumerate(unique_locations.iterrows()):
            if idx % 100 == 0:  # Log progress
                self.log(f"Processing SSRD correction for location {idx+1}/{total_locations}")
            
            lat, lon = location['latitude'], location['longitude']
            location_df = df_sorted[
                (df_sorted['latitude'] == lat) & (df_sorted['longitude'] == lon)
            ].copy()
            
            if len(location_df) == 0:
                continue
                
            location_df = location_df.sort_values('valid_time')
            
            # Extract hour
            location_df['hour'] = location_df['valid_time'].dt.hour
            
            # Calculate differences between consecutive SSRD values
            location_df['ssrd_diff'] = location_df['ssrd'].diff().fillna(0)
            
            # Apply correction logic
            location_df['ssrd_hourly'] = 0.0
            
            for i in range(len(location_df)):
                hour = location_df.iloc[i]['hour']
                ssrd_val = location_df.iloc[i]['ssrd']
                ssrd_diff = location_df.iloc[i]['ssrd_diff']
                
                if hour == 1:  # 01:00 UTC (first hour after midnight reset)
                    location_df.iloc[i, location_df.columns.get_loc('ssrd_hourly')] = ssrd_val
                elif ssrd_diff < 0 and hour != 1:  # Handle resets that aren't at 01:00
                    location_df.iloc[i, location_df.columns.get_loc('ssrd_hourly')] = 0
                else:
                    location_df.iloc[i, location_df.columns.get_loc('ssrd_hourly')] = ssrd_diff
            
            # Keep in J/m² - no conversion needed
            location_df['ssrd_corrected'] = location_df['ssrd_hourly']
            
            # Remove intermediate columns
            location_df = location_df.drop(columns=['hour', 'ssrd_diff', 'ssrd_hourly'])
            
            corrected_dfs.append(location_df)
        
        # Combine all corrected data
        result_df = pd.concat(corrected_dfs, ignore_index=True)
        
        # Ensure we maintain the original order
        result_df = result_df.sort_values(['valid_time', 'latitude', 'longitude'])
        
        self.log(f"SSRD correction complete. Original range: [{df['ssrd'].min():.2f}, {df['ssrd'].max():.2f}] J/m²")
        if 'ssrd_corrected' in result_df.columns:
            self.log(f"Corrected range: [{result_df['ssrd_corrected'].min():.2f}, {result_df['ssrd_corrected'].max():.2f}] J/m²")
        
        return result_df
    
    def process(self):
        """Run the complete processing workflow with checkpoint detection"""
        self.log(f"==================================================")
        self.log(f"Processing summer months for year {self.year}")
        self.log(f"==================================================")
        
        # Detect current step
        start_step = self.detect_current_step()
        
        if start_step >= 8:
            self.log("All processing steps already complete!")
            return True
        
        # Initialize variables for step dependencies
        grib_files = None
        csv_files = None
        analysis_results = None
        merged_files = None
        final_file = None
        
        # Step 1: Download or locate GRIB files
        if start_step <= 1:
            self.log("STEP 1: Download/locate GRIB files")
            grib_files = self.step1_get_grib_files()
            if not grib_files:
                self.log("Error: Failed to get GRIB files")
                return False
        else:
            # Load existing GRIB inventory
            grib_inventory_file = os.path.join(self.intermediate_dir, 'grib_inventory.json')
            if os.path.exists(grib_inventory_file):
                with open(grib_inventory_file, 'r') as f:
                    grib_data = json.load(f)
                    grib_files = {
                        'era5': {int(k): v for k, v in grib_data['era5'].items()},
                        'era5_land': {int(k): v for k, v in grib_data['era5_land'].items()}
                    }
                    self.log("Loaded existing GRIB file inventory")
        
        # Step 2: Convert GRIB files to CSV
        if start_step <= 2:
            self.log("STEP 2: Convert GRIB files to CSV")
            csv_files = self.step2_convert_grib_to_csv(grib_files)
            if not csv_files:
                self.log("Error: Failed to convert GRIB files to CSV")
                return False
        else:
            # Load existing CSV inventory
            csv_inventory_file = os.path.join(self.intermediate_dir, 'csv_inventory.json')
            if os.path.exists(csv_inventory_file):
                with open(csv_inventory_file, 'r') as f:
                    csv_files = json.load(f)
                    self.log("Loaded existing CSV file inventory")
        
        # Step 3: Analyze CSV files and save diagnostics
        if start_step <= 3:
            self.log("STEP 3: Analyze CSV files")
            analysis_results = self.step3_analyze_csv_files(csv_files)
            if not analysis_results:
                self.log("Error: Failed to analyze CSV files")
                return False
        else:
            self.log("STEP 3: CSV analysis already complete")
        
        # Step 4: Merge and interpolate data
        if start_step <= 4:
            self.log("STEP 4: Merge and interpolate data")
            merged_files = self.step4_merge_and_interpolate(csv_files['era5_rad'], csv_files['era5_land'])
            if not merged_files:
                self.log("Error: Failed to merge and interpolate data")
                return False
        else:
            # Load existing merged data
            merged_csv = os.path.join(self.intermediate_dir, "merged_interpolated.csv")
            if os.path.exists(merged_csv):
                merged_df = pd.read_csv(merged_csv)
                merged_df['valid_time'] = pd.to_datetime(merged_df['valid_time'])
                merged_files = {'merged_csv': merged_csv, 'merged_df': merged_df}
                self.log("Loaded existing merged data")
        
        # Step 5: Create final NetCDF file
        if start_step <= 5:
            self.log("STEP 5: Create final NetCDF file")
            final_file = self.step5_create_final_netcdf(merged_files)
            if not final_file:
                self.log("Error: Failed to create final NetCDF file")
                return False
        else:
            final_file = os.path.join(self.merged_dir, f"era5_summer_{self.year}.nc")
            self.log("Final NetCDF file already exists")
        
        # Step 6: Verify final file
        if start_step <= 6:
            self.log("STEP 6: Verify final file")
            self.step6_verify_final_file(final_file)
        
        # Step 7: Create SSRD and TP corrected version
        if start_step <= 7:
            self.log("STEP 7: Create SSRD and TP corrected version")
            corrected_file = self.step7_create_corrected_version(merged_files if merged_files else {'merged_csv': os.path.join(self.intermediate_dir, "merged_interpolated.csv")})
            if not corrected_file:
                self.log("Error: Failed to create corrected version")
                return False
        
        self.log(f"Successfully processed summer data for {self.year}")
        self.log(f"Final NetCDF file: {final_file}")
        self.log(f"SSRD and TP corrected files available in: {self.merged_corrected_dir}")
        
        return True
    
    def step7_create_corrected_version(self, merged_files):
        """Step 7: Create SSRD and TP corrected version of the data"""
        self.log("Creating SSRD and TP corrected version of the data...")
        
        # Load merged data if not already loaded
        merged_csv = merged_files.get('merged_csv')
        if not merged_csv or not os.path.exists(merged_csv):
            merged_csv = os.path.join(self.intermediate_dir, "merged_interpolated.csv")
            if not os.path.exists(merged_csv):
                self.log("Error: No merged data found for correction")
                return None
        
        self.log(f"Loading merged data from: {merged_csv}")
        df = pd.read_csv(merged_csv)
        df['valid_time'] = pd.to_datetime(df['valid_time'])
        
        # Apply SSRD correction
        df_corrected = self.correct_ssrd_joules(df)
        
        # Apply TP correction
        df_corrected = self.correct_tp(df_corrected)
        
        # Update ssrd and tp columns with corrected values
        if 'ssrd_corrected' in df_corrected.columns:
            df_corrected['ssrd'] = df_corrected['ssrd_corrected']
            df_corrected = df_corrected.drop(columns=['ssrd_corrected'])
        
        if 'tp_corrected' in df_corrected.columns:
            df_corrected['tp'] = df_corrected['tp_corrected']
            df_corrected = df_corrected.drop(columns=['tp_corrected'])
        
        # Save corrected CSV
        corrected_csv = os.path.join(self.merged_corrected_dir, f"era5_summer_{self.year}_corrected.csv")
        df_corrected.to_csv(corrected_csv, index=False)
        self.log(f"Saved corrected CSV: {corrected_csv}")
        
        # Create corrected NetCDF
        self.log("Creating corrected NetCDF file...")
        
        # Get unique coordinates and times
        unique_lons = sorted(df_corrected['longitude'].unique())
        unique_lats = sorted(df_corrected['latitude'].unique())
        unique_times = sorted(df_corrected['unix_time'].unique())
        
        # Create xarray Dataset
        ds = xr.Dataset(
            coords={
                'longitude': unique_lons,
                'latitude': unique_lats,
                'valid_time': unique_times
            }
        )
        
        # Add each variable in the correct order
        for var in self.variable_order:
            if var in df_corrected.columns:
                self.log(f"Adding corrected variable: {var}")
                
                # Initialize with NaN
                data_array = np.full((len(unique_lons), len(unique_lats), len(unique_times)), np.nan)
                
                # Create lookup dictionaries for faster indexing
                time_lookup = {t: i for i, t in enumerate(unique_times)}
                lat_lookup = {lat: i for i, lat in enumerate(unique_lats)}
                lon_lookup = {lon: i for i, lon in enumerate(unique_lons)}
                
                # Fill the array with data
                for _, row in df_corrected.iterrows():
                    lon_idx = lon_lookup.get(row['longitude'])
                    lat_idx = lat_lookup.get(row['latitude'])
                    time_idx = time_lookup.get(row['unix_time'])
                    
                    if lon_idx is not None and lat_idx is not None and time_idx is not None:
                        data_array[lon_idx, lat_idx, time_idx] = row[var]
                
                # Add to dataset
                ds[var] = xr.DataArray(
                    data_array,
                    dims=('longitude', 'latitude', 'valid_time'),
                    coords={
                        'longitude': unique_lons,
                        'latitude': unique_lats,
                        'valid_time': unique_times
                    }
                )
                
                # Add units
                if var == 't2m' or var == 'd2m':
                    ds[var].attrs['units'] = 'K'
                elif var == 'sp':
                    ds[var].attrs['units'] = 'Pa'
                elif var in ['u10', 'v10']:
                    ds[var].attrs['units'] = 'm s-1'
                elif var == 'tp':
                    ds[var].attrs['units'] = 'm'
                    ds[var].attrs['note'] = 'Corrected from cumulative m to hourly m'
                elif var == 'ssrd':
                    ds[var].attrs['units'] = 'J m-2'
                    ds[var].attrs['note'] = 'Corrected from cumulative J/m² to hourly J/m²'
                elif var == 'fdir':
                    ds[var].attrs['units'] = 'J m-2'
                elif var == 'tcc' or var == 'lsm':
                    ds[var].attrs['units'] = '(0 - 1)'
                elif var in ['msdwlwrf', 'msnlwrf']:
                    ds[var].attrs['units'] = 'W m-2'
        
        # Add metadata
        ds.attrs['title'] = f"ERA5 and ERA5-Land merged dataset for {self.year} (SSRD and tp corrected)"
        ds.attrs['source'] = "ERA5 and ERA5-Land processed for compatibility with NicheMapR micro_era5 function"
        ds.attrs['creation_date'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        ds.attrs['area'] = str(self.area)
        ds.attrs['processing_workflow'] = "Combined ERA5 radiation and ERA5-Land with interpolation and SSRD correction (J/m²)"
        ds.attrs['ssrd_correction'] = "SSRD converted from cumulative J/m² to hourly J/m²"
        ds.attrs['tp_correction'] = "tp converted from cumulative m to hourly m"
        
        # Set time units
        ds['valid_time'].attrs['units'] = 'seconds since 1970-01-01'
        
        # Save corrected NetCDF
        corrected_nc = os.path.join(self.merged_corrected_dir, f"era5_summer_{self.year}_corrected.nc")
        self.log(f"Saving corrected NetCDF: {corrected_nc}")
        
        ds.to_netcdf(
            corrected_nc,
            format='NETCDF3_CLASSIC',
            engine='netcdf4',
            encoding={
                'valid_time': {'dtype': 'int32'}
            }
        )
        

        
        return corrected_nc

    # Keep the original download and conversion methods (slightly modified for better logging)
    def download_era5_radiation(self, year, month):
        """Download ERA5 radiation data (from original script)"""
        month_str = f"{month:02d}"
        grib_file = os.path.join(self.era5_dir, f"era5_rad_{year}_{month_str}.grib")
        
        if os.path.exists(grib_file) and self.verify_grib_file(grib_file):
            self.log(f"ERA5 radiation data for {year}-{month_str} already exists: {grib_file}")
            return grib_file
        
        self.log(f"Downloading ERA5 radiation data for {year}-{month_str}...")
        
        if self.test_mode and month == 6:
            days = [f"{day:02d}" for day in range(1, min(6, calendar.monthrange(year, month)[1] + 1))]
            self.log(f"Test mode: Downloading only first {len(days)} days of June")
        else:
            days = [f"{day:02d}" for day in range(1, calendar.monthrange(year, month)[1] + 1)]
        
        north, west, south, east = self.area
        request_area = [north + 0.25, west - 0.25, south - 0.25, east + 0.25]
        
        request = {
            "product_type": "reanalysis",
            "variable": list(self.era5_variables.keys()),
            "year": str(year),
            "month": month_str,
            "day": days,
            "time": [f"{hour:02d}:00" for hour in range(24)],
            "area": request_area,
            "grid": [0.25, 0.25],
            "format": "grib"
        }
        
        try:
            self.log(f"Sending ERA5 radiation data request to CDS API...")
            downloaded_file = self.client.retrieve("reanalysis-era5-single-levels", request).download()
            
            # Handle ZIP files if needed (same logic as original)
            if downloaded_file.lower().endswith('.zip') or self.is_zip_file(downloaded_file):
                self.log(f"Downloaded file is a ZIP archive, extracting...")
                extract_dir = os.path.join(self.era5_dir, f"temp_extract_{year}_{month_str}")
                os.makedirs(extract_dir, exist_ok=True)
                
                with zipfile.ZipFile(downloaded_file, 'r') as zip_ref:
                    zip_ref.extractall(extract_dir)
                
                extracted_files = []
                for root, _, files in os.walk(extract_dir):
                    for file in files:
                        if file.lower().endswith(('.grib', '.grb', '.grib2', '.grb2')):
                            extracted_files.append(os.path.join(root, file))
                
                if extracted_files:
                    shutil.copy(extracted_files[0], grib_file)
                    shutil.rmtree(extract_dir)
                else:
                    self.log(f"No GRIB files found in ZIP archive")
                    return None
            else:
                if downloaded_file != grib_file:
                    shutil.move(downloaded_file, grib_file)
            
            if self.verify_grib_file(grib_file):
                return grib_file
            else:
                self.log(f"Downloaded file verification failed: {grib_file}")
                return None
                
        except Exception as e:
            self.log(f"Error downloading ERA5 radiation data: {e}")
            return None
    
    def download_era5_land(self, year, month):
        """Download ERA5-Land data (from original script)"""
        month_str = f"{month:02d}"
        grib_file = os.path.join(self.era5land_dir, f"era5_land_{year}_{month_str}.grib")
        
        if os.path.exists(grib_file) and self.verify_grib_file(grib_file):
            self.log(f"ERA5-Land data for {year}-{month_str} already exists: {grib_file}")
            return grib_file
        
        self.log(f"Downloading ERA5-Land data for {year}-{month_str}...")
        
        if self.test_mode and month == 6:
            days = [f"{day:02d}" for day in range(1, min(6, calendar.monthrange(year, month)[1] + 1))]
            self.log(f"Test mode: Downloading only first {len(days)} days of June")
        else:
            days = [f"{day:02d}" for day in range(1, calendar.monthrange(year, month)[1] + 1)]
        
        north, west, south, east = self.area
        request_area = [north + 0.1, west - 0.1, south - 0.1, east + 0.1]
        
        request = {
            "variable": list(self.era5land_variables.keys()),
            "year": str(year),
            "month": month_str,
            "day": days,
            "time": [f"{hour:02d}:00" for hour in range(24)],
            "area": request_area,
            "grid": [0.1, 0.1],
            "format": "grib"
        }
        
        try:
            self.log(f"Sending ERA5-Land data request to CDS API...")
            downloaded_file = self.client.retrieve("reanalysis-era5-land", request).download()
            
            # Handle ZIP files if needed (same logic as original)
            if downloaded_file.lower().endswith('.zip') or self.is_zip_file(downloaded_file):
                self.log(f"Downloaded file is a ZIP archive, extracting...")
                extract_dir = os.path.join(self.era5land_dir, f"temp_extract_{year}_{month_str}")
                os.makedirs(extract_dir, exist_ok=True)
                
                with zipfile.ZipFile(downloaded_file, 'r') as zip_ref:
                    zip_ref.extractall(extract_dir)
                
                extracted_files = []
                for root, _, files in os.walk(extract_dir):
                    for file in files:
                        if file.lower().endswith(('.grib', '.grb', '.grib2', '.grb2')):
                            extracted_files.append(os.path.join(root, file))
                
                if extracted_files:
                    shutil.copy(extracted_files[0], grib_file)
                    shutil.rmtree(extract_dir)
                else:
                    self.log(f"No GRIB files found in ZIP archive")
                    return None
            else:
                if downloaded_file != grib_file:
                    shutil.move(downloaded_file, grib_file)
            
            if self.verify_grib_file(grib_file):
                return grib_file
            else:
                self.log(f"Downloaded file verification failed: {grib_file}")
                return None
                
        except Exception as e:
            self.log(f"Error downloading ERA5-Land data: {e}")
            return None
    
    def is_zip_file(self, filepath):
        """Check if a file is a ZIP archive"""
        try:
            with open(filepath, 'rb') as f:
                magic = f.read(4)
                return magic.startswith(b'PK\x03\x04')
        except Exception:
            return False
    
    def verify_grib_file(self, grib_file):
        """Verify that a GRIB file is valid"""
        if not os.path.exists(grib_file) or os.path.getsize(grib_file) == 0:
            return False
        
        try:
            datasets = cfgrib.open_datasets(grib_file)
            return len(datasets) > 0
        except Exception:
            try:
                ds = xr.open_dataset(grib_file, engine='cfgrib')
                ds.close()
                return True
            except Exception:
                return False
    
    def filter_month_boundaries(self, df, month, year):
        """Filter dataframe to only include data from the specified month/year"""
        if 'valid_time' not in df.columns:
            return df
        
        df['valid_time'] = pd.to_datetime(df['valid_time'])
        
        start_date = pd.Timestamp(f"{year}-{month:02d}-01")
        if month == 12:
            end_date = pd.Timestamp(f"{year+1}-01-01") - pd.Timedelta(seconds=1)
        else:
            end_date = pd.Timestamp(f"{year}-{month+1:02d}-01") - pd.Timedelta(seconds=1)
        
        original_len = len(df)
        df_filtered = df[(df['valid_time'] >= start_date) & (df['valid_time'] <= end_date)]
        filtered_len = len(df_filtered)
        
        if original_len != filtered_len:
            self.log(f"Filtered month {month}/{year}: {filtered_len} rows kept, {original_len - filtered_len} rows removed")
        
        return df_filtered
    
    def convert_grib_to_csv(self, grib_file):
        """Convert a GRIB file to CSV with proper month boundary filtering"""
        self.log(f"Converting GRIB file to CSV: {grib_file}")
        
        try:
            filename = os.path.basename(grib_file)
            match = re.search(r'(\d{4})_(\d{2})', filename)
            if match:
                year = int(match.group(1))
                month = int(match.group(2))
                date_str = f"{match.group(2)}_{match.group(1)}"
            else:
                date_str = "unknown_date"
                self.log(f"Warning: Could not extract date from filename: {filename}")
                return []
            
            if "era5_rad" in filename:
                file_type = "era5_rad"
            elif "era5_land" in filename:
                file_type = "era5_land"
            else:
                file_type = "unknown"
            
            file_size_mb = os.path.getsize(grib_file) / (1024 * 1024)
            self.log(f"GRIB file size: {file_size_mb:.2f} MB")
            
            try:
                datasets = cfgrib.open_datasets(grib_file)
            except Exception as e:
                self.log(f"Error opening with cfgrib.open_datasets: {e}")
                try:
                    ds = xr.open_dataset(grib_file, engine='cfgrib')
                    datasets = [ds]
                except Exception as e2:
                    self.log(f"Error with xarray fallback: {e2}")
                    return []
            
            if not datasets:
                self.log(f"No datasets found in {grib_file}")
                return []
            
            csv_files = []
            for i, ds in enumerate(datasets):
                if len(datasets) > 1:
                    csv_file = os.path.join(self.csv_dir, f"{file_type}_{date_str}_part{i+1}.csv")
                else:
                    csv_file = os.path.join(self.csv_dir, f"{file_type}_{date_str}.csv")
                
                try:
                    df = ds.to_dataframe().reset_index()
                    
                    # Fix column names
                    if 'latitude' not in df.columns and 'lat' in df.columns:
                        df = df.rename(columns={'lat': 'latitude'})
                    if 'longitude' not in df.columns and 'lon' in df.columns:
                        df = df.rename(columns={'lon': 'longitude'})
                    
                    # Round coordinates
                    if 'latitude' in df.columns:
                        df['latitude'] = df['latitude'].round(2)
                    if 'longitude' in df.columns:
                        df['longitude'] = df['longitude'].round(2)
                    
                    # Create valid_time
                    if 'valid_time' not in df.columns:
                        if 'time' in df.columns and 'step' in df.columns:
                            try:
                                df['valid_time'] = pd.to_datetime(df['time']) + pd.to_timedelta(df['step'])
                            except:
                                df['valid_time'] = pd.to_datetime(df['time'])
                        elif 'time' in df.columns:
                            df['valid_time'] = pd.to_datetime(df['time'])
                    
                    # Filter month boundaries
                    df = self.filter_month_boundaries(df, month, year)
                    
                    if len(df) == 0:
                        self.log(f"Warning: No data left after filtering for {month}/{year}")
                        continue
                    
                    # Log time range
                    if 'valid_time' in df.columns:
                        time_min = df['valid_time'].min()
                        time_max = df['valid_time'].max()
                        self.log(f"Time range: {time_min} to {time_max}")
                    
                    df.to_csv(csv_file, index=False)
                    self.log(f"Created: {csv_file} (shape: {df.shape})")
                    csv_files.append(csv_file)
                    
                except Exception as e:
                    self.log(f"Error converting dataset to DataFrame: {e}")
                    continue
            
            return csv_files
            
        except Exception as e:
            self.log(f"Error converting GRIB to CSV: {e}")
            import traceback
            traceback.print_exc()
            return []


def parse_args():
    parser = argparse.ArgumentParser(description='Process ERA5 and ERA5-Land data with improved workflow')
    parser.add_argument('--years', type=int, nargs='+', help='Years to process (e.g. --years 2015 2016 2017)')
    parser.add_argument('--base-dir', type=str, default='era5_workflow_new', help='Base directory for all files')
    parser.add_argument('--test', action='store_true', help='Run in test mode (process only 5 days of June)')
    parser.add_argument('--force-overwrite', action='store_true', help='Force overwrite existing files and restart from beginning')
    return parser.parse_args()

def main():
    args = parse_args()
    
    if not args.years:
        years_input = input("Enter years to process (separated by spaces, e.g., '2015 2016 2017'): ")
        years = [int(year) for year in years_input.split()]
    else:
        years = args.years
    
    base_dir = args.base_dir
    os.makedirs(base_dir, exist_ok=True)
    
    for year in years:
        print(f"\n===== Processing year {year} =====")
        processor = ERA5CombinedProcessor(
            base_dir=base_dir, 
            year=year, 
            test_mode=args.test, 
            force_overwrite=args.force_overwrite
        )
        success = processor.process()
        
        if success:
            print(f"Successfully processed year {year}")
        else:
            print(f"Failed to process year {year}")
    
    print(f"\nAll years processed.")
    print(f"Original results are in: {os.path.join(base_dir, 'merged')}")
    print(f"SSRD-corrected results are in: {os.path.join(base_dir, 'merged_corrected')}")

if __name__ == "__main__":
    main()