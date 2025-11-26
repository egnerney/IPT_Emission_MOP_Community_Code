#!/usr/bin/env python3
"""
convert_npz_to_h5.py
====================
Convert Jovian Plasma Torus Model from .npz to HDF5 format

This script converts the compressed numpy .npz file to HDF5 format for
cross-platform compatibility with Python, IDL, Fortran, C++, and MATLAB.

The HDF5 format provides:
- Cross-language compatibility
- Built-in compression
- Metadata storage
- Self-documenting structure

Array Storage Convention:
- Arrays are stored in C-order (row-major, numpy default)
- IDL users will need to transpose arrays when reading
- Coordinate arrays: 1D vectors
- Data arrays: 3D with shape (nx, ny, nz) in C-order

Author: Edward (Eddie) G. Nerney Nov 2025
License: MIT
Version: 1.2 (Fixed array name mapping and metadata)
"""

import numpy as np
import h5py
import sys
from pathlib import Path


def convert_npz_to_h5(npz_file, h5_file, compression='gzip', compression_level=4):
    """
    Convert .npz plasma model to HDF5 format.
    
    Parameters
    ----------
    npz_file : str or Path
        Input .npz filename
    h5_file : str or Path
        Output .h5 filename
    compression : str, optional
        HDF5 compression filter ('gzip', 'lzf', or None)
        Default is 'gzip' for maximum compatibility
    compression_level : int, optional
        Compression level for gzip (0-9, higher = more compression)
        Default is 4 (balanced speed/size)
    
    Returns
    -------
    dict
        Dictionary of array names and shapes for verification
    """
    
    print("="*70)
    print("Converting Jovian Plasma Model: .npz → HDF5")
    print("="*70)
    print()
    
    # Load the .npz file
    print(f"Loading: {npz_file}")
    try:
        data = np.load(npz_file)
    except FileNotFoundError:
        print(f"Error: File not found: {npz_file}")
        print("\nPlease ensure the .npz file is in the current directory or provide full path.")
        return None
    
    print(f"Loaded {len(data.files)} arrays from .npz file")
    print()
    
    # Display contents
    print("Contents of .npz file:")
    print("-"*70)
    total_size_mb = 0
    array_info = {}
    
    for name in data.files:
        arr = data[name]
        size_mb = arr.nbytes / (1024**2)
        total_size_mb += size_mb
        array_info[name] = {
            'shape': arr.shape,
            'dtype': arr.dtype,
            'size_mb': size_mb
        }
        print(f"  {name:20s}  shape: {str(arr.shape):20s}  dtype: {str(arr.dtype):10s}  "
              f"size: {size_mb:8.2f} MB")
    
    print(f"\nTotal uncompressed size: {total_size_mb:.2f} MB")
    print()
    
    # Map array names to match expected names in code
    name_mapping = {
        'x_axis': 'x',
        'y_axis': 'y', 
        'z_axis': 'z',
        'nec': 'ne_c',
        'neh': 'ne_h',
        'Tec': 'Te_c',
        'Teh': 'Te_h',
        'Toph': 'To_ph',
        'Thp': 'T_hp'
    }
    
    # Create HDF5 file
    print(f"Creating HDF5 file: {h5_file}")
    print(f"Compression: {compression if compression else 'None'}")
    if compression == 'gzip':
        print(f"Compression level: {compression_level}")
    print()
    
    with h5py.File(h5_file, 'w') as hf:
        
        # Store global metadata
        hf.attrs['title'] = 'Jovian Plasma Torus Model - Io Plasma Torus'
        hf.attrs['description'] = 'Three-dimensional plasma model interpolated to regular grid'
        hf.attrs['source'] = 'Converted from jovian_plasma_interpolated_381x381x231.npz'
        hf.attrs['units_system'] = 'Jupiter radii (R_J), cm^-3, eV'
        hf.attrs['coordinate_system'] = 'Jovian System III (right-handed)'
        hf.attrs['array_order'] = 'C-order (row-major, numpy default)'
        hf.attrs['note_for_idl_users'] = 'Arrays must be transposed when reading into IDL'
        hf.attrs['note_for_fortran_users'] = 'Arrays must be transposed when reading into Fortran'
        hf.attrs['creation_date'] = '2025-11-07'
        hf.attrs['version'] = '1.2'
        
        # Separate coordinate arrays from data arrays
        coord_group = hf.create_group('coordinates')
        data_group = hf.create_group('data')
        
        # Identify coordinate arrays (1D) vs data arrays (3D)
        coord_arrays = {}
        data_arrays = {}
        
        for name in data.files:
            arr = data[name]
            if arr.ndim == 1:
                coord_arrays[name] = arr
            else:
                data_arrays[name] = arr
        
        # Store coordinate arrays (no compression needed for 1D)
        print("Storing coordinate arrays:")
        print("-"*70)
        for name, arr in coord_arrays.items():
            output_name = name_mapping.get(name, name)  # Apply mapping
            coord_group.create_dataset(output_name, data=arr, dtype=arr.dtype)
            
            # Add metadata for coordinate arrays
            if output_name == 'x':
                coord_group[output_name].attrs['description'] = 'X coordinate (equatorial plane)'
                coord_group[output_name].attrs['units'] = 'R_J (Jupiter radii)'
                coord_group[output_name].attrs['axis'] = 'X (anti-sunward positive)'
            elif output_name == 'y':
                coord_group[output_name].attrs['description'] = 'Y coordinate (equatorial plane)'
                coord_group[output_name].attrs['units'] = 'R_J (Jupiter radii)'
                coord_group[output_name].attrs['axis'] = 'Y (dusk positive)'
            elif output_name == 'z':
                coord_group[output_name].attrs['description'] = 'Z coordinate (rotation axis)'
                coord_group[output_name].attrs['units'] = 'R_J (Jupiter radii)'
                coord_group[output_name].attrs['axis'] = 'Z (north positive)'
            
            print(f"  {name:20s}  [{len(arr)}]  → /coordinates/{output_name}")
        
        print()
        
        # Store data arrays with compression
        print("Storing data arrays with compression:")
        print("-"*70)
        
        compression_opts = None
        if compression == 'gzip':
            compression_opts = compression_level
        
        for name, arr in data_arrays.items():
            output_name = name_mapping.get(name, name)  # Apply mapping
            # Store in C-order (numpy default)
            data_group.create_dataset(
                output_name, 
                data=arr, 
                dtype=arr.dtype,
                compression=compression,
                compression_opts=compression_opts,
                shuffle=True  # Improves compression
            )
            
            # Add metadata for data arrays
            if output_name == 'ne_c':
                data_group[output_name].attrs['description'] = 'Electron number density (cold component)'
                data_group[output_name].attrs['units'] = 'cm^-3'
                data_group[output_name].attrs['component'] = 'cold'
            elif output_name == 'Te_c':
                data_group[output_name].attrs['description'] = 'Electron temperature (cold component)'
                data_group[output_name].attrs['units'] = 'eV'
                data_group[output_name].attrs['component'] = 'cold'
            elif output_name == 'ne_h':
                data_group[output_name].attrs['description'] = 'Electron number density (hot component)'
                data_group[output_name].attrs['units'] = 'cm^-3'
                data_group[output_name].attrs['component'] = 'hot'
            elif output_name == 'Te_h':
                data_group[output_name].attrs['description'] = 'Electron temperature (hot component)'
                data_group[output_name].attrs['units'] = 'eV'
                data_group[output_name].attrs['component'] = 'hot'
            elif output_name == 'nsp':
                data_group[output_name].attrs['description'] = 'S+ (sulfur single ionized) number density'
                data_group[output_name].attrs['units'] = 'cm^-3'
                data_group[output_name].attrs['ion'] = 'S+'
            elif output_name == 'ns2p':
                data_group[output_name].attrs['description'] = 'S++ (sulfur double ionized) number density'
                data_group[output_name].attrs['units'] = 'cm^-3'
                data_group[output_name].attrs['ion'] = 'S++'
            elif output_name == 'ns3p':
                data_group[output_name].attrs['description'] = 'S+++ (sulfur triple ionized) number density'
                data_group[output_name].attrs['units'] = 'cm^-3'
                data_group[output_name].attrs['ion'] = 'S+++'
            elif output_name == 'nop':
                data_group[output_name].attrs['description'] = 'O+ (oxygen single ionized) number density'
                data_group[output_name].attrs['units'] = 'cm^-3'
                data_group[output_name].attrs['ion'] = 'O+'
            elif output_name == 'no2p':
                data_group[output_name].attrs['description'] = 'O++ (oxygen double ionized) number density'
                data_group[output_name].attrs['units'] = 'cm^-3'
                data_group[output_name].attrs['ion'] = 'O++'
            elif output_name == 'noph':
                data_group[output_name].attrs['description'] = 'O+ hot component number density'
                data_group[output_name].attrs['units'] = 'cm^-3'
                data_group[output_name].attrs['ion'] = 'O+ (hot)'
            elif output_name == 'nhp':
                data_group[output_name].attrs['description'] = 'H+ (hydrogen ionized) number density'
                data_group[output_name].attrs['units'] = 'cm^-3'
                data_group[output_name].attrs['ion'] = 'H+'
            elif output_name == 'nnap':
                data_group[output_name].attrs['description'] = 'Na+ (sodium ionized) number density'
                data_group[output_name].attrs['units'] = 'cm^-3'
                data_group[output_name].attrs['ion'] = 'Na+'
            elif output_name == 'Ti':
                data_group[output_name].attrs['description'] = 'Ion temperature'
                data_group[output_name].attrs['units'] = 'eV'
            elif output_name == 'To_ph':
                data_group[output_name].attrs['description'] = 'O+ hot component temperature'
                data_group[output_name].attrs['units'] = 'eV'
            elif output_name == 'T_hp':
                data_group[output_name].attrs['description'] = 'H+ temperature'
                data_group[output_name].attrs['units'] = 'eV'
            
            data_group[output_name].attrs['shape_description'] = '(nx, ny, nz) in C-order'
            data_group[output_name].attrs['dimensions'] = f"{arr.shape[0]} x {arr.shape[1]} x {arr.shape[2]}"
            
            original_size = arr.nbytes / (1024**2)
            print(f"  {name:20s}  {str(arr.shape):20s}  {original_size:8.2f} MB  "
                  f"→ /data/{output_name}")
        
        print()
    
    # Get file sizes
    npz_size = Path(npz_file).stat().st_size / (1024**2)
    h5_size = Path(h5_file).stat().st_size / (1024**2)
    
    print("="*70)
    print("Conversion Complete")
    print("="*70)
    print(f"Original .npz file:     {npz_size:8.2f} MB")
    print(f"New .h5 file:           {h5_size:8.2f} MB")
    print(f"Compression ratio:      {npz_size/h5_size:.2f}x")
    print(f"Size reduction:         {(1 - h5_size/npz_size)*100:.1f}%")
    print()
    print("File structure:")
    print("  /coordinates/         (1D coordinate arrays)")
    print("  /data/                (3D plasma parameter arrays)")
    print()
    print(f"Output file: {h5_file}")
    print()
    
    return array_info


def verify_h5_file(h5_file):
    """
    Verify the contents of the HDF5 file.
    
    Parameters
    ----------
    h5_file : str or Path
        HDF5 filename to verify
    """
    print("="*70)
    print("Verifying HDF5 File Structure")
    print("="*70)
    print()
    
    with h5py.File(h5_file, 'r') as hf:
        # Print global attributes
        print("Global Attributes:")
        print("-"*70)
        for key, value in hf.attrs.items():
            print(f"  {key:25s}: {value}")
        print()
        
        # Print coordinate arrays
        print("Coordinate Arrays (/coordinates/):")
        print("-"*70)
        for name in hf['coordinates'].keys():
            dset = hf['coordinates'][name]
            print(f"  {name:20s}  shape: {dset.shape}  dtype: {dset.dtype}")
            if dset.attrs:
                for attr_name, attr_value in dset.attrs.items():
                    print(f"    {attr_name}: {attr_value}")
        print()
        
        # Print data arrays
        print("Data Arrays (/data/):")
        print("-"*70)
        for name in hf['data'].keys():
            dset = hf['data'][name]
            print(f"  {name:20s}  shape: {dset.shape}  dtype: {dset.dtype}")
            print(f"    compression: {dset.compression}")
            if 'description' in dset.attrs:
                print(f"    description: {dset.attrs['description']}")
            if 'units' in dset.attrs:
                print(f"    units: {dset.attrs['units']}")
        print()


def main():
    """Main conversion function."""
    
    # Check command line arguments
    if len(sys.argv) > 1:
        npz_file = sys.argv[1]
    else:
        npz_file = 'jovian_plasma_interpolated_381x381x231.npz'
    
    if len(sys.argv) > 2:
        h5_file = sys.argv[2]
    else:
        # Generate output filename
        h5_file = Path(npz_file).stem + '.h5'
    
    # Convert
    array_info = convert_npz_to_h5(
        npz_file, 
        h5_file, 
        compression='gzip', 
        compression_level=4
    )
    
    if array_info is None:
        return 1
    
    # Verify
    print()
    verify_h5_file(h5_file)
    
    print("="*70)
    print("Conversion successful!")
    print("="*70)
    print()
    print("Usage in Python:")
    print("  import h5py")
    print(f"  with h5py.File('{h5_file}', 'r') as f:")
    print("      x = f['coordinates/x'][:]")
    print("      ne_c = f['data/ne_c'][:]")
    print()
    print("Usage in IDL:")
    print(f"  file_id = H5F_OPEN('{h5_file}')")
    print("  x = H5D_READ(H5D_OPEN(file_id, '/coordinates/x'))")
    print("  ne_c = H5D_READ(H5D_OPEN(file_id, '/data/ne_c'))")
    print("  ne_c = TRANSPOSE(ne_c)  ; Transpose for IDL row-major order")
    print("  H5F_CLOSE, file_id")
    print()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())