#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
convert_hdf5_strings.py

HDF5 Variable-Length String Conversion Script for Fortran Compatibility
========================================================================

This script converts HDF5 files containing variable-length (VL) string datasets
(Python h5py default format) to use fixed-length ASCII strings that can be read
by Fortran's standard HDF5 bindings.

PROBLEM:
Python h5py writes strings as variable-length (VL) strings by default, stored
using HDF5's hvl_t structure. Fortran's standard HDF5 Fortran bindings cannot
read these VL strings directly, producing "unable to convert between src and
dest datatype" errors.

SOLUTION:
This script rewrites the 'species' dataset in the emission table HDF5 files
to use fixed-length ASCII strings (dtype 'S4' for 4-byte fixed strings),
which Fortran can read directly using H5T_FORTRAN_S1 memory type.

USAGE:
    python convert_hdf5_strings.py

The script will process all emission table files in ../Emiss_tables/ relative
to its location, creating Fortran-compatible versions.

AUTHOR: Edward (Eddie) G. Nerney
INSTITUTION: Laboratory for Atmospheric and Space Physics, CU Boulder
LICENSE: Open source for academic and research use
DATE: November 2025
"""

import h5py
import numpy as np
from pathlib import Path
import shutil
import sys


def check_string_dataset_type(filename, dataset_name):
    """
    Check the string type of a dataset in an HDF5 file.
    
    Parameters
    ----------
    filename : str or Path
        Path to HDF5 file
    dataset_name : str
        Name of the dataset to check
        
    Returns
    -------
    str
        Description of the string type
    """
    with h5py.File(filename, 'r') as f:
        if dataset_name not in f:
            return "Dataset not found"
        
        dset = f[dataset_name]
        dtype = dset.dtype
        
        if dtype.kind == 'O':
            return "Variable-length string (object dtype) - NOT Fortran compatible"
        elif dtype.kind == 'S':
            return f"Fixed-length string ({dtype.itemsize} bytes) - Fortran compatible"
        elif h5py.check_string_dtype(dtype):
            info = h5py.check_string_dtype(dtype)
            if info.length is None:
                return "Variable-length string (vlen) - NOT Fortran compatible"
            else:
                return f"Fixed-length string ({info.length} bytes) - Fortran compatible"
        else:
            return f"Non-string type: {dtype}"


def convert_species_to_fixed_length(input_file, output_file=None, max_length=4):
    """
    Convert the 'species' dataset from variable-length to fixed-length strings.
    
    Parameters
    ----------
    input_file : str or Path
        Path to input HDF5 file with VL string 'species' dataset
    output_file : str or Path, optional
        Path to output file. If None, overwrites input file (after backup).
    max_length : int, optional
        Maximum string length for fixed-length format (default: 4)
        
    Returns
    -------
    bool
        True if conversion successful, False otherwise
    """
    input_file = Path(input_file)
    
    if not input_file.exists():
        print(f"ERROR: Input file not found: {input_file}")
        return False
    
    # Determine output file
    if output_file is None:
        # Create backup and overwrite original
        backup_file = input_file.with_suffix('.h5.bak')
        if not backup_file.exists():
            print(f"  Creating backup: {backup_file.name}")
            shutil.copy2(input_file, backup_file)
        output_file = input_file
        in_place = True
    else:
        output_file = Path(output_file)
        in_place = False
    
    print(f"Processing: {input_file.name}")
    
    # Read original file
    try:
        with h5py.File(input_file, 'r') as f:
            if 'species' not in f:
                print("  No 'species' dataset found - skipping")
                return True
            
            # Check current type
            dtype_info = check_string_dataset_type(input_file, 'species')
            print(f"  Original species dtype: {dtype_info}")
            
            if "Fortran compatible" in dtype_info and "NOT" not in dtype_info:
                print("  Already Fortran compatible - no conversion needed")
                return True
            
            # Read all data
            species_vl = f['species'][:]
            
            # Convert to Python strings, handling various formats
            species_strings = []
            for s in species_vl:
                if isinstance(s, bytes):
                    species_strings.append(s.decode('utf-8', errors='replace'))
                elif isinstance(s, str):
                    species_strings.append(s)
                else:
                    species_strings.append(str(s))
            
            # Find actual max length needed
            actual_max = max(len(s) for s in species_strings)
            fixed_len = max(max_length, actual_max)
            print(f"  Max species string length: {actual_max}, using fixed length: {fixed_len}")
            
            # Convert to fixed-length ASCII
            species_fixed = np.array(species_strings, dtype=f'S{fixed_len}')
            
            # Collect all other datasets
            all_data = {}
            all_attrs = dict(f.attrs)
            
            def collect_data(name, obj):
                if isinstance(obj, h5py.Dataset):
                    if name != 'species':
                        all_data[name] = {
                            'data': obj[:],
                            'dtype': obj.dtype,
                            'attrs': dict(obj.attrs)
                        }
                    else:
                        all_data[name] = {
                            'data': species_fixed,
                            'dtype': species_fixed.dtype,
                            'attrs': dict(obj.attrs)
                        }
            
            f.visititems(collect_data)
    
    except Exception as e:
        print(f"  ERROR reading file: {e}")
        return False
    
    # Write new file
    try:
        # If in-place, use temporary file first
        if in_place:
            temp_file = input_file.with_suffix('.h5.tmp')
            write_file = temp_file
        else:
            write_file = output_file
        
        with h5py.File(write_file, 'w') as f:
            # Write file attributes
            for key, val in all_attrs.items():
                f.attrs[key] = val
            
            # Write all datasets
            for name, info in all_data.items():
                # Handle grouped datasets (with '/' in name)
                if '/' in name:
                    group_name = name.rsplit('/', 1)[0]
                    if group_name not in f:
                        f.create_group(group_name)
                
                # Create dataset with explicit dtype
                dset = f.create_dataset(name, data=info['data'], dtype=info['dtype'])
                
                # Copy attributes
                for key, val in info['attrs'].items():
                    dset.attrs[key] = val
        
        # If in-place, replace original with temp
        if in_place:
            temp_file.replace(input_file)
        
        # Verify conversion
        new_dtype_info = check_string_dataset_type(output_file, 'species')
        print(f"  New species dtype: {new_dtype_info}")
        
        if "Fortran compatible" in new_dtype_info and "NOT" not in new_dtype_info:
            print("  Conversion successful!")
            return True
        else:
            print("  WARNING: Conversion may not have worked as expected")
            return False
    
    except Exception as e:
        print(f"  ERROR writing file: {e}")
        return False


def main():
    """
    Main function to convert all emission table files.
    """
    print("="*70)
    print("HDF5 Variable-Length String Conversion for Fortran Compatibility")
    print("="*70)
    print()
    
    # Determine paths
    script_dir = Path(__file__).parent
    
    # Try different relative paths to find Emiss_tables
    possible_paths = [
        script_dir / ".." / "Emiss_tables",
        script_dir / "Emiss_tables",
        script_dir.parent.parent / "Emiss_tables",
        Path("../Emiss_tables"),
        Path("Emiss_tables"),
    ]
    
    emiss_dir = None
    for path in possible_paths:
        if path.exists():
            emiss_dir = path.resolve()
            break
    
    if emiss_dir is None:
        print("ERROR: Could not find Emiss_tables directory")
        print("Please ensure the directory structure is:")
        print("  IPT_Emission_MOP_Community_Code/")
        print("    LOS_Integration/")
        print("      convert_hdf5_strings.py  (this script)")
        print("    Emiss_tables/")
        print("      CHIANTI_11.0.2_emiss_tables_*.h5")
        return 1
    
    print(f"Found Emiss_tables directory: {emiss_dir}")
    print()
    
    # Find emission table files
    h5_files = list(emiss_dir.glob("CHIANTI*emiss_tables*.h5"))
    
    if not h5_files:
        print("No emission table HDF5 files found in Emiss_tables directory")
        return 1
    
    print(f"Found {len(h5_files)} emission table file(s):")
    for f in h5_files:
        print(f"  - {f.name}")
    print()
    
    # Convert each file
    success_count = 0
    for h5_file in h5_files:
        if convert_species_to_fixed_length(h5_file):
            success_count += 1
        print()
    
    print("="*70)
    print(f"Conversion complete: {success_count}/{len(h5_files)} files processed successfully")
    print("="*70)
    
    if success_count == len(h5_files):
        print("\nAll files are now Fortran-compatible!")
        print("You can now compile and run the Fortran code.")
        return 0
    else:
        print("\nSome files may need manual attention.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
