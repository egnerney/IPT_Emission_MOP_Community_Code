#!/usr/bin/env python3
"""
convert_idl_emiss_tables_sav_files_to_h5.py

Convert CHIANTI IDL .sav emission table files to HDF5 format.
Tested on Python 3.12, NumPy 2.3.4, SciPy 1.16.3

Usage:
    python convert_sav_to_h5.py

This script will:
1. Read IDL .sav files from ../Emiss_tables/
2. Convert float64 → float32 (50% size reduction)
3. Save as HDF5 files with proper naming (lowercase n, uppercase T)
4. Handle IDL column-major → NumPy row-major ordering
"""

import sys
from pathlib import Path
import numpy as np
import h5py
from scipy.io import readsav
from datetime import datetime


def inspect_structure(struct_array, name=""):
    """
    Recursively inspect IDL structure arrays.
    
    Parameters
    ----------
    struct_array : np.ndarray of structs
        IDL structure array
    name : str
        Name for printing
        
    Returns
    -------
    dict
        Flattened dictionary of arrays from structure
    """
    result = {}
    
    if not hasattr(struct_array, 'dtype') or struct_array.dtype.names is None:
        return result
    
    # Get the first element if it's an array
    if struct_array.shape:
        struct = struct_array.flat[0]
    else:
        struct = struct_array
    
    # Extract all fields
    for field_name in struct_array.dtype.names:
        field_data = struct_array[field_name]
        
        # If it's a single element array, extract it
        if field_data.shape == (1,):
            field_data = field_data[0]
        
        # If it's nested structures, recurse
        if hasattr(field_data, 'dtype') and field_data.dtype.names:
            nested = inspect_structure(field_data, f"{name}.{field_name}")
            result.update(nested)
        else:
            key = f"{name}.{field_name}" if name else field_name
            result[key] = field_data
    
    return result


def convert_single_maxwellian(sav_file, output_file):
    """
    Convert single Maxwellian .sav file to HDF5.
    
    Parameters
    ----------
    sav_file : str or Path
        Path to input .sav file
    output_file : str or Path
        Path to output .h5 file
    """
    print("\n" + "="*60)
    print("Converting Single Maxwellian Tables")
    print("="*60)
    print(f"Reading: {sav_file}")
    
    # Read IDL save file
    data = readsav(str(sav_file), python_dict=True, verbose=False)
    
    print(f"\nVariables in .sav file:")
    for key in sorted(data.keys()):
        value = data[key]
        if isinstance(value, np.ndarray):
            print(f"  {key:<20} {str(value.shape):<20} {value.dtype}")
        else:
            print(f"  {key:<20} {type(value).__name__}")
    
    # Extract temperature and density grids
    if 'temp_arr' in data:
        T = data['temp_arr'].astype(np.float32)
    elif 'temperature' in data:
        T = data['temperature'].astype(np.float32)
    else:
        raise ValueError("Cannot find temperature array (tried 'temp_arr', 'temperature')")
    
    if 'dens_arr' in data:
        n = data['dens_arr'].astype(np.float32)
    elif 'density' in data:
        n = data['density'].astype(np.float32)
    else:
        raise ValueError("Cannot find density array (tried 'dens_arr', 'density')")
    
    print(f"\n✓ Temperature grid: {T.shape}, range [{T.min():.2e}, {T.max():.2e}] eV")
    print(f"✓ Density grid: {n.shape}, range [{n.min():.2e}, {n.max():.2e}] cm^-3")
    
    # Extract emissivity data organized by species from YPTSI structure
    # IDL saves as (n_T, n_n, n_lines) but NumPy reads as (n_lines, n_n, n_T)
    # Need to transpose to (n_T, n_n, n_lines)
    
    species_names = ['S', 'SP', 'S2P', 'S3P', 'S4P', 'O', 'OP', 'O2P', 'NAP']
    emiss_arrays = []
    wave_arrays = []
    species_labels = []
    
    if 'yptsi' in data and 'xwavi' in data:
        print("\n✓ Found YPTSI and XWAVI structures")
        print("  Extracting species-organized emission data...")
        
        yptsi = data['yptsi'][0]  # Structure array with 1 element
        xwavi = data['xwavi'][0]
        
        for species_name in species_names:
            # Get field name (IDL structures use uppercase)
            field_name = species_name.upper()
            
            # Check if field exists (try both cases)
            if hasattr(yptsi, field_name):
                emiss_data = getattr(yptsi, field_name)
            elif hasattr(yptsi, species_name):
                emiss_data = getattr(yptsi, species_name)
            else:
                print(f"  ⚠ Skipping {species_name} (not found)")
                continue
            
            # Get wavelength data
            if hasattr(xwavi, field_name):
                wave_data = getattr(xwavi, field_name)
            elif hasattr(xwavi, species_name):
                wave_data = getattr(xwavi, species_name)
            else:
                print(f"  ⚠ Skipping {species_name} (no wavelength data)")
                continue
            
            # IDL saves as (n_T, n_n, n_lines)
            # NumPy reads as (n_lines, n_n, n_T) due to column vs row major
            # Transpose to (n_T, n_n, n_lines)
            if emiss_data.ndim == 3:
                emiss_data = emiss_data.transpose(2, 1, 0)
                print(f"  ✓ {species_name:5s}: {emiss_data.shape[2]:5d} lines, "
                      f"transposed to {emiss_data.shape}")
            else:
                print(f"  ⚠ {species_name}: unexpected shape {emiss_data.shape}")
                continue
            
            # Store
            emiss_arrays.append(emiss_data.astype(np.float32))
            wave_arrays.append(wave_data.astype(np.float32))
            species_labels.extend([species_name] * len(wave_data))
        
        # Concatenate all species along wavelength dimension
        emissivity = np.concatenate(emiss_arrays, axis=2)  # Concat along line dimension
        wavelength = np.concatenate(wave_arrays)
        species_array = np.array(species_labels, dtype='U10')
        
        print(f"\n✓ Combined emission data:")
        print(f"  Total lines: {emissivity.shape[2]}")
        print(f"  Final shape: {emissivity.shape} (T, n, lines)")
        print(f"  Species breakdown:")
        for name in species_names:
            count = sum(1 for s in species_labels if s == name)
            if count > 0:
                print(f"    {name:5s}: {count:5d} lines")
    
    else:
        raise ValueError("Could not find YPTSI and XWAVI structures in .sav file")
    
    print(f"\n✓ Emissivity array: {emissivity.shape}")
    
    # Create HDF5 file
    print(f"\nWriting HDF5 file: {output_file}")
    with h5py.File(output_file, 'w') as f:
        # Write grids
        f.create_dataset('T', data=T, compression='gzip', compression_opts=4)
        f.create_dataset('n', data=n, compression='gzip', compression_opts=4)
        
        # Write emissivity
        f.create_dataset('emissivity', data=emissivity, compression='gzip', compression_opts=4)
        
        # Write wavelength
        f.create_dataset('wavelength', data=wavelength, compression='gzip', compression_opts=4)
        print(f"  ✓ Wavelength: {len(wavelength)} values")
        
        # Write species
        # Convert to bytes for h5py compatibility
        species_bytes = [s.encode('utf-8') for s in species_array]
        dt = h5py.string_dtype(encoding='utf-8')
        f.create_dataset('species', data=species_bytes, dtype=dt)
        print(f"  ✓ Species: {len(species_array)} entries")
        
        # Add attributes
        f.attrs['description'] = 'CHIANTI 11.0.2 single Maxwellian emission tables'
        f.attrs['version'] = '11.0.2'
        f.attrs['creation_date'] = datetime.now().isoformat()
        f.attrs['grid_size'] = f'{len(T)}x{len(n)}'
        f.attrs['data_type'] = 'float32'
        f.attrs['original_file'] = str(sav_file)
        
        # Dataset attributes
        f['T'].attrs['description'] = 'Temperature grid (eV)'
        f['T'].attrs['units'] = 'eV'
        f['n'].attrs['description'] = 'Density grid (cm^-3)'
        f['n'].attrs['units'] = 'cm^-3'
        f['emissivity'].attrs['description'] = 'Emissivity tables (erg cm^3 s^-1)'
        f['emissivity'].attrs['units'] = 'erg cm^3 s^-1'
        f['wavelength'].attrs['description'] = 'Wavelength (Angstrom)'
        f['wavelength'].attrs['units'] = 'Angstrom'
    
    # Report file sizes
    sav_size = Path(sav_file).stat().st_size / (1024**2)
    h5_size = Path(output_file).stat().st_size / (1024**2)
    
    print(f"\n✓ Conversion complete!")
    print(f"  Original .sav: {sav_size:.1f} MB")
    print(f"  New .h5: {h5_size:.1f} MB")
    print(f"  Savings: {100*(1-h5_size/sav_size):.1f}%")


def convert_double_maxwellian(sav_file, output_file):
    """
    Convert double Maxwellian .sav file to HDF5.
    
    Parameters
    ----------
    sav_file : str or Path
        Path to input .sav file
    output_file : str or Path
        Path to output .h5 file
    """
    print("\n" + "="*60)
    print("Converting Double Maxwellian Tables")
    print("="*60)
    print(f"Reading: {sav_file}")
    
    # Read IDL save file
    data = readsav(str(sav_file), python_dict=True, verbose=False)
    
    print(f"\nVariables in .sav file:")
    for key in sorted(data.keys()):
        value = data[key]
        if isinstance(value, np.ndarray):
            print(f"  {key:<20} {str(value.shape):<20} {value.dtype}")
        else:
            print(f"  {key:<20} {type(value).__name__}")
    
    # Extract grids based on actual variable names
    # From your file: ne_arr, tec_arr, teh_arr, feh_arr
    
    if 'tec_arr' in data:
        T_cold = data['tec_arr'].astype(np.float32)
    else:
        raise ValueError("Could not find tec_arr (cold electron temperature)")
    
    if 'teh_arr' in data:
        T_hot = data['teh_arr'].astype(np.float32)
    else:
        raise ValueError("Could not find teh_arr (hot electron temperature)")
    
    if 'ne_arr' in data:
        n = data['ne_arr'].astype(np.float32)  # Total electron density
    elif 'ne_total_arr' in data:
        n = data['ne_total_arr'].astype(np.float32)
    else:
        raise ValueError("Could not find ne_arr (total electron density)")
    
    if 'feh_arr' in data:
        feh = data['feh_arr'].astype(np.float32)  # Hot electron fraction
    else:
        raise ValueError("Could not find feh_arr (hot electron fraction)")
    
    print(f"\n✓ T_cold: {T_cold.shape}, range [{T_cold.min():.2e}, {T_cold.max():.2e}] eV")
    print(f"✓ T_hot: {T_hot.shape}, range [{T_hot.min():.2e}, {T_hot.max():.2e}] eV")
    print(f"✓ n_total: {n.shape}, range [{n.min():.2e}, {n.max():.2e}] cm^-3")
    print(f"✓ feh: {feh.shape}, range [{feh.min():.2e}, {feh.max():.2e}]")
    
    # Extract emissivity data organized by species from YPTSI structure
    # IDL saves as (n_Tc, n_Th, n_n, n_feh, n_lines)
    # NumPy reads as (n_lines, n_feh, n_n, n_Th, n_Tc) due to column vs row major
    # Need to transpose to (n_Tc, n_Th, n_n, n_feh, n_lines)
    
    species_names = ['S', 'SP', 'S2P', 'S3P', 'S4P', 'O', 'OP', 'O2P', 'NAP']
    emiss_arrays = []
    wave_arrays = []
    species_labels = []
    
    if 'yptsi' in data and 'xwavi' in data:
        print("\n✓ Found YPTSI and XWAVI structures")
        print("  Extracting species-organized emission data...")
        
        yptsi = data['yptsi'][0]  # Structure array with 1 element
        xwavi = data['xwavi'][0]
        
        for species_name in species_names:
            # Get field name (IDL structures use uppercase)
            field_name = species_name.upper()
            
            # Check if field exists (try both cases)
            if hasattr(yptsi, field_name):
                emiss_data = getattr(yptsi, field_name)
            elif hasattr(yptsi, species_name):
                emiss_data = getattr(yptsi, species_name)
            else:
                print(f"  ⚠ Skipping {species_name} (not found)")
                continue
            
            # Get wavelength data
            if hasattr(xwavi, field_name):
                wave_data = getattr(xwavi, field_name)
            elif hasattr(xwavi, species_name):
                wave_data = getattr(xwavi, species_name)
            else:
                print(f"  ⚠ Skipping {species_name} (no wavelength data)")
                continue
            
            # IDL saves as (n_Tc, n_Th, n_n, n_feh, n_lines)
            # NumPy reads as (n_lines, n_feh, n_n, n_Th, n_Tc) due to dimension reversal
            # Transpose to (n_Tc, n_Th, n_n, n_feh, n_lines)
            if emiss_data.ndim == 5:
                emiss_data = emiss_data.transpose(4, 3, 2, 1, 0)
                print(f"  ✓ {species_name:5s}: {emiss_data.shape[4]:5d} lines, "
                      f"transposed to {emiss_data.shape}")
            else:
                print(f"  ⚠ {species_name}: unexpected shape {emiss_data.shape}")
                continue
            
            # Store
            emiss_arrays.append(emiss_data.astype(np.float32))
            wave_arrays.append(wave_data.astype(np.float32))
            species_labels.extend([species_name] * len(wave_data))
        
        # Concatenate all species along wavelength dimension
        emissivity = np.concatenate(emiss_arrays, axis=4)  # Concat along line dimension
        wavelength = np.concatenate(wave_arrays)
        species_array = np.array(species_labels, dtype='U10')
        
        print(f"\n✓ Combined emission data:")
        print(f"  Total lines: {emissivity.shape[4]}")
        print(f"  Final shape: {emissivity.shape} (T_cold, T_hot, n, feh, lines)")
        print(f"  Species breakdown:")
        for name in species_names:
            count = sum(1 for s in species_labels if s == name)
            if count > 0:
                print(f"    {name:5s}: {count:5d} lines")
    
    else:
        raise ValueError("Could not find YPTSI and XWAVI structures in .sav file")
    
    print(f"\n✓ Emissivity array: {emissivity.shape}")
    
    # Create HDF5 file
    print(f"\nWriting HDF5 file: {output_file}")
    with h5py.File(output_file, 'w') as f:
        # Write grids
        f.create_dataset('T_cold', data=T_cold, compression='gzip', compression_opts=4)
        f.create_dataset('T_hot', data=T_hot, compression='gzip', compression_opts=4)
        f.create_dataset('n', data=n, compression='gzip', compression_opts=4)
        f.create_dataset('feh', data=feh, compression='gzip', compression_opts=4)
        
        # Write emissivity
        f.create_dataset('emissivity', data=emissivity, compression='gzip', compression_opts=4)
        
        # Write wavelength
        f.create_dataset('wavelength', data=wavelength, compression='gzip', compression_opts=4)
        print(f"  ✓ Wavelength: {len(wavelength)} values")
        
        # Write species
        # Convert to bytes for h5py compatibility
        species_bytes = [s.encode('utf-8') for s in species_array]
        dt = h5py.string_dtype(encoding='utf-8')
        f.create_dataset('species', data=species_bytes, dtype=dt)
        print(f"  ✓ Species: {len(species_array)} entries")
        
        # Add attributes
        f.attrs['description'] = 'CHIANTI 11.0.2 double Maxwellian emission tables'
        f.attrs['version'] = '11.0.2'
        f.attrs['creation_date'] = datetime.now().isoformat()
        f.attrs['grid_size'] = f'{len(T_cold)}x{len(T_hot)}x{len(n)}x{len(feh)}'
        f.attrs['data_type'] = 'float32'
        f.attrs['original_file'] = str(sav_file)
        
        # Dataset attributes
        f['T_cold'].attrs['description'] = 'Cold electron temperature grid (eV)'
        f['T_cold'].attrs['units'] = 'eV'
        f['T_hot'].attrs['description'] = 'Hot electron temperature grid (eV)'
        f['T_hot'].attrs['units'] = 'eV'
        f['n'].attrs['description'] = 'Total electron density grid (cm^-3)'
        f['n'].attrs['units'] = 'cm^-3'
        f['feh'].attrs['description'] = 'Hot electron fraction (neh/(nec+neh))'
        f['feh'].attrs['units'] = 'dimensionless'
        f['emissivity'].attrs['description'] = 'Emissivity tables (erg cm^3 s^-1)'
        f['emissivity'].attrs['units'] = 'erg cm^3 s^-1'
        f['wavelength'].attrs['description'] = 'Wavelength (Angstrom)'
        f['wavelength'].attrs['units'] = 'Angstrom'
    
    # Report file sizes
    sav_size = Path(sav_file).stat().st_size / (1024**2)
    h5_size = Path(output_file).stat().st_size / (1024**2)
    
    print(f"\n✓ Conversion complete!")
    print(f"  Original .sav: {sav_size:.1f} MB")
    print(f"  New .h5: {h5_size:.1f} MB")
    print(f"  Savings: {100*(1-h5_size/sav_size):.1f}%")


def main():
    """Main conversion routine."""
    
    print("="*60)
    print("CHIANTI Emission Tables: .sav → .h5 Converter")
    print("="*60)
    print(f"Python: {sys.version}")
    print(f"NumPy: {np.__version__}")
    print(f"h5py: {h5py.__version__}")
    
    # Define file paths
    base_dir = Path(__file__).parent
    emiss_dir = base_dir / '..' / 'Emiss_tables'
    
    # Single Maxwellian file
    single_sav = emiss_dir / 'CHIANTI_11.0.2_emiss_arrays_all_species_all_wavelengths_50x50_logspaced.sav'
    single_h5 = emiss_dir / 'CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5'
    
    # Double Maxwellian file
    double_sav = emiss_dir / 'CHIANTI_11.0.2_emiss_arrays_all_species_all_wavelengths_15x10x20x10_hote_logspaced.sav'
    double_h5 = emiss_dir / 'CHIANTI_11.0.2_emiss_tables_double_maxwellian_15x10x20x10.h5'
    
    # Convert single Maxwellian
    if single_sav.exists():
        try:
            convert_single_maxwellian(single_sav, single_h5)
        except Exception as e:
            print(f"\n✗ Error converting single Maxwellian: {e}")
            import traceback
            traceback.print_exc()
    else:
        print(f"\n⚠ Single Maxwellian file not found: {single_sav}")
    
    # Convert double Maxwellian
    if double_sav.exists():
        try:
            convert_double_maxwellian(double_sav, double_h5)
        except Exception as e:
            print(f"\n✗ Error converting double Maxwellian: {e}")
            import traceback
            traceback.print_exc()
    else:
        print(f"\n⚠ Double Maxwellian file not found: {double_sav}")
    
    print("\n" + "="*60)
    print("Conversion Summary")
    print("="*60)
    
    if single_h5.exists():
        print(f"✓ Single Maxwellian: {single_h5.name}")
    if double_h5.exists():
        print(f"✓ Double Maxwellian: {double_h5.name}")
    


if __name__ == '__main__':
    main()