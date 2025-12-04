classdef IPT_emiss_MOP_community_code < handle
    % IPT_emiss_MOP_community_code
    %
    % Io Plasma Torus (IPT) UV/Optical Emission Line-of-Sight Integration Module
    % ==========================================================================
    %
    % This class provides comprehensive functionality for calculating UV and optical
    % emission spectra from the Io Plasma Torus (IPT) using ray tracing through a 3D
    % plasma model with CHIANTI atomic emission tables. Supports both single and double
    % Maxwellian electron distributions with proper species-by-species interpolation
    % and vectorized computations.
    %
    % PHYSICAL MODEL:
    % - Plasma torus extends from ~5-10 R_J in cylindrical radius
    % - Peak emission near 6 R_J at centrifugal equator
    % - Scale height ~0.5-1 R_J
    % - Based on IPT Isotropic Model interpolated to rectilinear grid
    % - Per-ion photon emission rates from CHIANTI 11.0.2 atomic database
    % - Line-of-sight integration using Simpson's rule
    % - Optically thin plasma approximation (valid for IPT)
    %
    % COMPUTATIONAL APPROACH:
    % - Ray tracing through 3D plasma model with trilinear interpolation
    % - Bilinear interpolation of single Maxwellian emission tables (2D: Te, ne)
    % - Quadrilinear interpolation of double Maxwellian emission tables (4D: Tec, Teh, ne, feh)
    % - Vectorized emission rate interpolation for all lines simultaneously
    % - Per-ion photon emission rates multiplied by ion density and 1e-6 Rayleigh factor
    % - Simpson's rule integration over line of sight
    % - Analytic ERF-based Gaussian convolution for instrument response
    %
    % COORDINATE SYSTEMS AND UNITS:
    % - Positions: Jupiter radii [R_J]
    % - Temperatures: electron volts [eV]
    % - Densities: particles per cubic centimeter [cm^-3]
    % - Wavelengths: Angstroms [A]
    % - Emission rates: photons per second per ion [photons s^-1 ion^-1]
    % - Brightnesses: Rayleighs [R], where 1 R = 10^6 photons s^-1 cm^-2 (4*pi sr)^-1
    %
    % APPLICABLE WAVELENGTH RANGES:
    % - UV instruments: 550-2100 A (JUICE-UVS, Europa-UVS, HST/STIS)
    % - Optical instruments: 3000-10000 A (ground-based telescopes, HST optical)
    %
    % REFERENCES:
    % - CHIANTI database: Dere et al. 1997; Del Zanna et al. 2020; Dufresne et al. 2024
    % - IPT observations: Steffl et al. 2004a,b; Thomas et al. 2004; Bagenal & Delamere 2011
    % - Emission modeling: Nerney et al. 2017, 2020, 2022, 2025a, 2025b
    % - Electron distributions: Meyer-Vernet & Moncuquet 1989; Moncuquet et al. 2002
    %
    % AUTHOR: Edward (Eddie) G. Nerney
    % INSTITUTION: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
    % LICENSE: Open source for academic and research use
    % VERSION: 1.0
    % DATE: November 2025
    %
    % CHIANTI ACKNOWLEDGMENT:
    % CHIANTI is a collaborative project involving George Mason University, the University
    % of Michigan (USA), University of Cambridge (UK), and NASA Goddard Space Flight Center (USA).
    
    properties (Constant)
        % Physical constants
        R_J = 71492.0;              % Jupiter radius in km
        R_J_CM = 7.1492e9;          % Jupiter radius in cm
        RAYLEIGH_FACTOR = 1e-6 * 7.1492e9;  % Conversion factor to Rayleighs
    end
    
    properties
        % Coordinate axes [R_J]
        x_axis
        y_axis
        z_axis
        
        % Electron density fields [cm^-3]
        nec         % Cold electron density
        neh         % Hot electron density
        ne_total    % Total electron density
        feh         % Hot electron fraction
        
        % Ion density fields [cm^-3]
        nsp         % S+ density
        ns2p        % S++ density
        ns3p        % S+++ density
        nop         % O+ density
        no2p        % O++ density
        
        % Electron temperature fields [eV]
        Tec         % Cold electron temperature
        Teh         % Hot electron temperature
        
        % 3D Interpolators for plasma parameters
        interp_nec
        interp_neh
        interp_ne_total
        interp_feh
        interp_nsp
        interp_ns2p
        interp_ns3p
        interp_nop
        interp_no2p
        interp_Tec
        interp_Teh
        
        % Single Maxwellian emission table parameters
        temp_arr        % Temperature grid [eV]
        dens_arr        % Density grid [cm^-3]
        log_temp        % Log10 of temperature grid
        log_dens        % Log10 of density grid
        
        % Double Maxwellian emission table parameters
        tec_arr         % Cold temperature grid [eV]
        teh_arr         % Hot temperature grid [eV]
        ne_arr          % Total density grid [cm^-3]
        feh_arr         % Hot fraction grid
        log_tec         % Log10 of cold temperature grid
        log_teh         % Log10 of hot temperature grid
        log_ne          % Log10 of density grid
        log_feh         % Log10 of hot fraction grid
        
        % Emission data organized by species (containers.Map)
        wavelengths_single      % Single Maxwellian wavelengths by species
        emissivities_single     % Single Maxwellian emission rates by species
        wavelengths_double      % Double Maxwellian wavelengths by species
        emissivities_double     % Double Maxwellian emission rates by species
        
        % State flags
        single_maxwellian_loaded = false
        double_maxwellian_loaded = false
        
        % Species display names mapping
        ion_names
        
        % Default species list for calculations
        default_species
        
        % Ion interpolator map for convenient access
        ion_interp_map
    end
    
    methods
        function obj = IPT_emiss_MOP_community_code(plasma_file, emission_file_single, emission_file_double)
            % Constructor - Initialize the raytracer with plasma model and emission tables.
            %
            % Parameters
            % ----------
            % plasma_file : char, optional
            %     Path to the plasma model HDF5 file. If not provided, uses default location.
            % emission_file_single : char, optional
            %     Path to single Maxwellian emission tables HDF5 file.
            % emission_file_double : char, optional
            %     Path to double Maxwellian emission tables HDF5 file.
            
            fprintf('Initializing Jovian UV/Optical Emission Raytracer...\n');
            
            % Species display names for output
            obj.ion_names = containers.Map(...
                {'SP', 'S2P', 'S3P', 'S4P', 'OP', 'O2P'}, ...
                {'S II', 'S III', 'S IV', 'S V', 'O II', 'O III'});
            
            % Default species for calculations (excludes trace species)
            obj.default_species = {'SP', 'S2P', 'S3P', 'OP', 'O2P'};
            
            % Initialize containers.Map for emission data storage
            obj.wavelengths_single = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.emissivities_single = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.wavelengths_double = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.emissivities_double = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.ion_interp_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
            
            % Set default paths if not provided
            if nargin < 1 || isempty(plasma_file)
                plasma_file = obj.get_default_plasma_path();
            end
            if nargin < 2 || isempty(emission_file_single)
                emission_file_single = obj.get_default_emission_path_single();
            end
            if nargin < 3 || isempty(emission_file_double)
                emission_file_double = obj.get_default_emission_path_double();
            end
            
            % Load all data files
            obj.load_plasma_model(plasma_file);
            obj.load_emission_tables_single(emission_file_single);
            obj.load_emission_tables_double(emission_file_double);
            
            fprintf('Initialization complete.\n');
        end
        
        function plasma_file = get_default_plasma_path(obj)
            % Get the default path to the plasma model file.
            module_dir = fileparts(mfilename('fullpath'));
            plasma_dir = fullfile(fileparts(fileparts(module_dir)), '3D_Torus_Model');
            plasma_file = fullfile(plasma_dir, 'jovian_plasma_interpolated_381x381x231.h5');
            
            if ~isfile(plasma_file)
                error('Plasma model not found at expected location: %s', plasma_file);
            end
        end
        
        function emiss_file = get_default_emission_path_single(obj)
            % Get the default path to the single Maxwellian emission tables file.
            module_dir = fileparts(mfilename('fullpath'));
            emiss_dir = fullfile(fileparts(fileparts(module_dir)), 'Emiss_tables');
            emiss_file = fullfile(emiss_dir, 'CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5');
            
            if ~isfile(emiss_file)
                error('Single Maxwellian emission tables not found: %s', emiss_file);
            end
        end
        
        function emiss_file = get_default_emission_path_double(obj)
            % Get the default path to the double Maxwellian emission tables file.
            module_dir = fileparts(mfilename('fullpath'));
            emiss_dir = fullfile(fileparts(fileparts(module_dir)), 'Emiss_tables');
            emiss_file = fullfile(emiss_dir, 'CHIANTI_11.0.2_emiss_tables_double_maxwellian_16x8x12x10.h5');
            
            if ~isfile(emiss_file)
                warning('Double Maxwellian emission tables not found: %s', emiss_file);
                emiss_file = '';
            end
        end
        
        function load_plasma_model(obj, filename)
            % Load 3D plasma model from HDF5 file.
            %
            % Parameters
            % ----------
            % filename : char
            %     Path to HDF5 file containing the 3D plasma model
            %
            % Notes
            % -----
            % The HDF5 file contains coordinate arrays and 3D fields for electron
            % densities, ion densities, and temperatures on a rectilinear grid.
            % Arrays are permuted from HDF5 storage order to MATLAB convention.
            
            fprintf('Loading plasma model from %s...\n', filename);
            
            % Load coordinate axes
            obj.x_axis = double(h5read(filename, '/coordinates/x'));
            obj.y_axis = double(h5read(filename, '/coordinates/y'));
            obj.z_axis = double(h5read(filename, '/coordinates/z'));
            obj.x_axis = obj.x_axis(:);
            obj.y_axis = obj.y_axis(:);
            obj.z_axis = obj.z_axis(:);
            
            % Load 3D plasma fields
            % HDF5/Python stores as (nx, ny, nz) but MATLAB reads as (nz, ny, nx)
            % Permute [3,2,1] converts back to (nx, ny, nz) for griddedInterpolant
            obj.nec = permute(double(h5read(filename, '/data/ne_c')), [3, 2, 1]);
            obj.neh = permute(double(h5read(filename, '/data/ne_h')), [3, 2, 1]);
            obj.nsp = permute(double(h5read(filename, '/data/nsp')), [3, 2, 1]);
            obj.ns2p = permute(double(h5read(filename, '/data/ns2p')), [3, 2, 1]);
            obj.ns3p = permute(double(h5read(filename, '/data/ns3p')), [3, 2, 1]);
            obj.nop = permute(double(h5read(filename, '/data/nop')), [3, 2, 1]);
            obj.no2p = permute(double(h5read(filename, '/data/no2p')), [3, 2, 1]);
            obj.Tec = permute(double(h5read(filename, '/data/Te_c')), [3, 2, 1]);
            obj.Teh = permute(double(h5read(filename, '/data/Te_h')), [3, 2, 1]);
            
            % Calculate derived quantities
            obj.ne_total = obj.nec + obj.neh;
            obj.feh = zeros(size(obj.ne_total));
            valid_ne = obj.ne_total > 0;
            obj.feh(valid_ne) = obj.neh(valid_ne) ./ obj.ne_total(valid_ne);
            
            % Clean up invalid values (NaN, negative, non-finite -> 0)
            obj.nec(~isfinite(obj.nec) | obj.nec < 0) = 0;
            obj.neh(~isfinite(obj.neh) | obj.neh < 0) = 0;
            obj.ne_total(~isfinite(obj.ne_total) | obj.ne_total < 0) = 0;
            obj.feh(~isfinite(obj.feh) | obj.feh < 0) = 0;
            obj.nsp(~isfinite(obj.nsp) | obj.nsp < 0) = 0;
            obj.ns2p(~isfinite(obj.ns2p) | obj.ns2p < 0) = 0;
            obj.ns3p(~isfinite(obj.ns3p) | obj.ns3p < 0) = 0;
            obj.nop(~isfinite(obj.nop) | obj.nop < 0) = 0;
            obj.no2p(~isfinite(obj.no2p) | obj.no2p < 0) = 0;
            obj.Tec(~isfinite(obj.Tec) | obj.Tec < 0) = 0;
            obj.Teh(~isfinite(obj.Teh) | obj.Teh < 0) = 0;
            
            % Create 3D interpolators
            fprintf('Creating 3D interpolators for plasma parameters...\n');
            obj.create_plasma_interpolators();
            
            fprintf('Plasma model loaded: grid shape [%d, %d, %d]\n', ...
                    size(obj.nec, 1), size(obj.nec, 2), size(obj.nec, 3));
            fprintf('  X range: %.1f to %.1f R_J\n', min(obj.x_axis), max(obj.x_axis));
            fprintf('  Y range: %.1f to %.1f R_J\n', min(obj.y_axis), max(obj.y_axis));
            fprintf('  Z range: %.1f to %.1f R_J\n', min(obj.z_axis), max(obj.z_axis));
        end
        
        function create_plasma_interpolators(obj)
            % Create 3D interpolators for all plasma parameters.
            %
            % Uses MATLAB's griddedInterpolant for trilinear interpolation
            % on non-uniform rectilinear grids. Points outside grid return NaN.
            
            % Create interpolators with 'linear' method and 'none' extrapolation
            obj.interp_nec = griddedInterpolant({obj.x_axis, obj.y_axis, obj.z_axis}, ...
                                                 obj.nec, 'linear', 'none');
            obj.interp_neh = griddedInterpolant({obj.x_axis, obj.y_axis, obj.z_axis}, ...
                                                 obj.neh, 'linear', 'none');
            obj.interp_ne_total = griddedInterpolant({obj.x_axis, obj.y_axis, obj.z_axis}, ...
                                                      obj.ne_total, 'linear', 'none');
            obj.interp_feh = griddedInterpolant({obj.x_axis, obj.y_axis, obj.z_axis}, ...
                                                 obj.feh, 'linear', 'none');
            obj.interp_nsp = griddedInterpolant({obj.x_axis, obj.y_axis, obj.z_axis}, ...
                                                 obj.nsp, 'linear', 'none');
            obj.interp_ns2p = griddedInterpolant({obj.x_axis, obj.y_axis, obj.z_axis}, ...
                                                  obj.ns2p, 'linear', 'none');
            obj.interp_ns3p = griddedInterpolant({obj.x_axis, obj.y_axis, obj.z_axis}, ...
                                                  obj.ns3p, 'linear', 'none');
            obj.interp_nop = griddedInterpolant({obj.x_axis, obj.y_axis, obj.z_axis}, ...
                                                 obj.nop, 'linear', 'none');
            obj.interp_no2p = griddedInterpolant({obj.x_axis, obj.y_axis, obj.z_axis}, ...
                                                  obj.no2p, 'linear', 'none');
            obj.interp_Tec = griddedInterpolant({obj.x_axis, obj.y_axis, obj.z_axis}, ...
                                                 obj.Tec, 'linear', 'none');
            obj.interp_Teh = griddedInterpolant({obj.x_axis, obj.y_axis, obj.z_axis}, ...
                                                 obj.Teh, 'linear', 'none');
            
            % Store ion interpolator mapping for convenient species-based access
            obj.ion_interp_map('SP') = obj.interp_nsp;
            obj.ion_interp_map('S2P') = obj.interp_ns2p;
            obj.ion_interp_map('S3P') = obj.interp_ns3p;
            obj.ion_interp_map('OP') = obj.interp_nop;
            obj.ion_interp_map('O2P') = obj.interp_no2p;
        end
        
        function load_emission_tables_single(obj, filename)
            % Load single Maxwellian emission tables from HDF5 file.
            %
            % Parameters
            % ----------
            % filename : char
            %     Path to HDF5 file containing single Maxwellian emission tables
            %
            % Table Structure
            % ---------------
            % The HDF5 file contains:
            % - T : temperature grid [eV], shape (n_T,)
            % - n : density grid [cm^-3], shape (n_n,)
            % - emiss : per-ion photon emission rates [photons s^-1 ion^-1]
            % - wavelength : emission line wavelengths [Angstroms]
            % - species : species labels for each line
            
            fprintf('Loading single Maxwellian emission tables from %s...\n', filename);
            
            % Load parameter grids
            obj.temp_arr = double(h5read(filename, '/T'));
            obj.dens_arr = double(h5read(filename, '/n'));
            obj.temp_arr = obj.temp_arr(:);
            obj.dens_arr = obj.dens_arr(:);
            
            % Load emission array
            % Python writes shape: (n_T, n_n, n_lines)
            % MATLAB h5read returns: (n_lines, n_n, n_T)
            emissivity_all = double(h5read(filename, '/emiss'));
            
            % Load wavelengths
            wavelength_all = double(h5read(filename, '/wavelength'));
            wavelength_all = wavelength_all(:);
            
            % Load and clean species labels (handle HDF5 string encoding)
            species_raw = h5read(filename, '/species');
            species_all = obj.clean_hdf5_strings(species_raw);
            
            % Create log-space grids for interpolation
            obj.log_temp = log10(obj.temp_arr);
            obj.log_dens = log10(obj.dens_arr);
            
            % Get unique species preserving order of first occurrence
            unique_species = {};
            for i = 1:length(species_all)
                s = species_all{i};
                if ~any(strcmp(unique_species, s))
                    unique_species{end+1} = s;
                end
            end
            
            % Organize data by species for efficient access
            for k = 1:length(unique_species)
                species_key = unique_species{k};
                indices = find(strcmp(species_all, species_key));
                
                obj.wavelengths_single(species_key) = wavelength_all(indices);
                
                % Extract emission data for this species
                % From MATLAB order (n_lines, n_n, n_T) extract subset then
                % permute to (n_lines_species, n_T, n_n) for interpolation
                emiss_species = emissivity_all(indices, :, :);
                obj.emissivities_single(species_key) = permute(emiss_species, [1, 3, 2]);
            end
            
            obj.single_maxwellian_loaded = true;
            
            fprintf('Single Maxwellian tables loaded successfully:\n');
            fprintf('  Temperature range: %.3f - %.1f eV\n', min(obj.temp_arr), max(obj.temp_arr));
            fprintf('  Density range: %.1f - %.0f cm^-3\n', min(obj.dens_arr), max(obj.dens_arr));
            fprintf('  Grid size: %d x %d\n', length(obj.temp_arr), length(obj.dens_arr));
            fprintf('  Species: %s\n', strjoin(sort(keys(obj.wavelengths_single)), ', '));
        end
        
        function load_emission_tables_double(obj, filename)
            % Load double Maxwellian emission tables from HDF5 file.
            %
            % Parameters
            % ----------
            % filename : char
            %     Path to HDF5 file containing double Maxwellian emission tables
            %
            % Physical Model
            % --------------
            % The double Maxwellian electron distribution function is:
            % f(v) = (1 - feh) * f_Maxwell(v, Tec) + feh * f_Maxwell(v, Teh)
            %
            % High-excitation lines are strongly enhanced by hot electrons.
            
            if isempty(filename)
                obj.double_maxwellian_loaded = false;
                return;
            end
            
            fprintf('Loading double Maxwellian emission tables from %s...\n', filename);
            
            % Load parameter grids
            obj.tec_arr = double(h5read(filename, '/T_cold'));
            obj.teh_arr = double(h5read(filename, '/T_hot'));
            obj.ne_arr = double(h5read(filename, '/n'));
            obj.feh_arr = double(h5read(filename, '/feh'));
            obj.tec_arr = obj.tec_arr(:);
            obj.teh_arr = obj.teh_arr(:);
            obj.ne_arr = obj.ne_arr(:);
            obj.feh_arr = obj.feh_arr(:);
            
            % Load emission array
            % Python writes shape: (n_Tc, n_Th, n_n, n_feh, n_lines)
            % MATLAB h5read returns: (n_lines, n_feh, n_n, n_Th, n_Tc)
            emissivity_all = double(h5read(filename, '/emiss'));
            
            % Load wavelengths
            wavelength_all = double(h5read(filename, '/wavelength'));
            wavelength_all = wavelength_all(:);
            
            % Load and clean species labels
            species_raw = h5read(filename, '/species');
            species_all = obj.clean_hdf5_strings(species_raw);
            
            % Create log-space grids for interpolation
            obj.log_tec = log10(obj.tec_arr);
            obj.log_teh = log10(obj.teh_arr);
            obj.log_ne = log10(obj.ne_arr);
            obj.log_feh = log10(obj.feh_arr);
            
            % Get unique species preserving order
            unique_species = {};
            for i = 1:length(species_all)
                s = species_all{i};
                if ~any(strcmp(unique_species, s))
                    unique_species{end+1} = s;
                end
            end
            
            % Organize data by species
            for k = 1:length(unique_species)
                species_key = unique_species{k};
                indices = find(strcmp(species_all, species_key));
                
                obj.wavelengths_double(species_key) = wavelength_all(indices);
                
                % Extract and permute emission data
                % From (n_lines, n_feh, n_n, n_Th, n_Tc)
                % To   (n_lines, n_Tc, n_Th, n_n, n_feh) for interpolation
                emiss_species = emissivity_all(indices, :, :, :, :);
                obj.emissivities_double(species_key) = permute(emiss_species, [1, 5, 4, 3, 2]);
            end
            
            obj.double_maxwellian_loaded = true;
            
            fprintf('Double Maxwellian tables loaded successfully:\n');
            fprintf('  Core temperature range: %.3f - %.1f eV\n', min(obj.tec_arr), max(obj.tec_arr));
            fprintf('  Hot temperature range: %.1f - %.0f eV\n', min(obj.teh_arr), max(obj.teh_arr));
            fprintf('  Density range: %.1f - %.0f cm^-3\n', min(obj.ne_arr), max(obj.ne_arr));
            fprintf('  Hot fraction range: %.6f - %.4f\n', min(obj.feh_arr), max(obj.feh_arr));
            fprintf('  Grid size: %d x %d x %d x %d\n', length(obj.tec_arr), length(obj.teh_arr), ...
                    length(obj.ne_arr), length(obj.feh_arr));
            fprintf('  Species: %s\n', strjoin(sort(keys(obj.wavelengths_double)), ', '));
        end
        
        function species_clean = clean_hdf5_strings(~, species_raw)
            % Clean species strings from HDF5 file.
            %
            % HDF5 string data can come in various formats depending on how it
            % was written and how MATLAB reads it. This function handles:
            % - 2D char arrays (fixed-width strings)
            % - Cell arrays of char vectors
            % - String arrays
            % - Byte arrays with null terminators
            %
            % Parameters
            % ----------
            % species_raw : various
            %     Raw species data from h5read
            %
            % Returns
            % -------
            % species_clean : cell array of char
            %     Cleaned species labels
            
            if ischar(species_raw)
                % 2D char array - each row is one string
                n = size(species_raw, 1);
                species_clean = cell(n, 1);
                for i = 1:n
                    s = species_raw(i, :);
                    % Remove null characters (char(0)) which terminate HDF5 strings
                    s = s(s ~= char(0));
                    % Remove any whitespace
                    s = strtrim(s);
                    species_clean{i} = s;
                end
            elseif iscell(species_raw)
                % Cell array of various types
                n = numel(species_raw);
                species_clean = cell(n, 1);
                for i = 1:n
                    s = species_raw{i};
                    if isnumeric(s)
                        % Byte array - convert to char
                        s = char(s(:)');
                    elseif isstring(s)
                        s = char(s);
                    end
                    % Remove nulls and whitespace
                    s = s(s ~= char(0));
                    s = strtrim(s);
                    species_clean{i} = s;
                end
            elseif isstring(species_raw)
                % String array
                species_clean = cellstr(species_raw);
                for i = 1:length(species_clean)
                    s = species_clean{i};
                    s = s(s ~= char(0));
                    species_clean{i} = strtrim(s);
                end
            else
                % Try to convert to char
                try
                    species_clean = cellstr(char(species_raw));
                    for i = 1:length(species_clean)
                        s = species_clean{i};
                        s = s(s ~= char(0));
                        species_clean{i} = strtrim(s);
                    end
                catch
                    error('Unable to parse species array of type: %s', class(species_raw));
                end
            end
            
            species_clean = species_clean(:);
        end
        
        function [s_values, positions] = trace_ray(obj, start_pos, direction, ds, max_distance)
            % Trace a ray through the plasma model.
            %
            % Parameters
            % ----------
            % start_pos : array, shape (3,)
            %     Starting position [x, y, z] in Jupiter radii
            % direction : array, shape (3,)
            %     Direction vector (will be normalized)
            % ds : float, optional
            %     Step size along ray in R_J (default: 0.1)
            % max_distance : float, optional
            %     Maximum distance to trace in R_J (default: 40.0)
            %
            % Returns
            % -------
            % s_values : array, shape (n_points,)
            %     Distances along ray in R_J
            % positions : array, shape (n_points, 3)
            %     Cartesian positions [x, y, z] along ray in R_J
            
            if nargin < 4 || isempty(ds)
                ds = 0.1;
            end
            if nargin < 5 || isempty(max_distance)
                max_distance = 40.0;
            end
            
            % Normalize direction vector
            direction = double(direction(:)');
            direction = direction / norm(direction);
            
            % Create ray points with uniform spacing
            n_points = floor(max_distance / ds) + 1;
            s_values = linspace(0.0, max_distance, n_points)';
            
            % Calculate positions along ray using broadcasting
            start_pos = double(start_pos(:)');
            positions = start_pos + s_values * direction;
        end
        
        function [wavelengths, emission_rates] = interpolate_emission_rates_single_vectorized(obj, ...
                Te_arr, ne_arr, species_key, min_wav, max_wav)
            % Interpolate per-ion photon emission rates for single Maxwellian.
            %
            % Uses vectorized bilinear interpolation in log10(T)-log10(n) space.
            %
            % Parameters
            % ----------
            % Te_arr : array, shape (n_los,)
            %     Electron temperatures along line of sight [eV]
            % ne_arr : array, shape (n_los,)
            %     Electron densities along line of sight [cm^-3]
            % species_key : char
            %     Species identifier ('SP', 'S2P', 'S3P', 'OP', 'O2P')
            % min_wav : float, optional
            %     Minimum wavelength to include [Angstroms]
            % max_wav : float, optional
            %     Maximum wavelength to include [Angstroms]
            %
            % Returns
            % -------
            % wavelengths : array, shape (n_lines,)
            %     Emission line wavelengths [Angstroms]
            % emission_rates : array, shape (n_los, n_lines)
            %     Per-ion photon emission rates [photons s^-1 ion^-1]
            
            if nargin < 5, min_wav = 550.0; end
            if nargin < 6, max_wav = 2100.0; end
            
            Te_arr = Te_arr(:);
            ne_arr = ne_arr(:);
            
            % Check if species exists
            if ~isKey(obj.wavelengths_single, species_key)
                wavelengths = [];
                emission_rates = zeros(length(Te_arr), 0);
                return;
            end
            
            % Get wavelengths and filter to range
            wav_all = obj.wavelengths_single(species_key);
            wav_mask = (wav_all >= min_wav) & (wav_all <= max_wav);
            wavelengths = wav_all(wav_mask);
            
            if isempty(wavelengths)
                wavelengths = [];
                emission_rates = zeros(length(Te_arr), 0);
                return;
            end
            
            % Get emission table: shape (n_lines, n_T, n_n)
            emiss_table_all = obj.emissivities_single(species_key);
            emiss_table = emiss_table_all(wav_mask, :, :);
            
            n_los = length(Te_arr);
            n_lines = length(wavelengths);
            emission_rates = zeros(n_los, n_lines);
            
            % Check for valid points
            valid = (Te_arr > 0) & (ne_arr > 0) & isfinite(Te_arr) & isfinite(ne_arr);
            in_bounds = valid & ...
                       (Te_arr >= obj.temp_arr(1)) & (Te_arr <= obj.temp_arr(end)) & ...
                       (ne_arr >= obj.dens_arr(1)) & (ne_arr <= obj.dens_arr(end));
            
            if ~any(in_bounds)
                return;
            end
            
            % Convert valid points to log space
            log_T = log10(Te_arr(in_bounds));
            log_n = log10(ne_arr(in_bounds));
            n_valid = sum(in_bounds);
            
            % Find bracketing indices using binary search
            i_T = obj.searchsorted_clip(obj.log_temp, log_T, 1, length(obj.log_temp) - 1);
            i_n = obj.searchsorted_clip(obj.log_dens, log_n, 1, length(obj.log_dens) - 1);
            
            i_T0 = i_T;
            i_T1 = i_T + 1;
            i_n0 = i_n;
            i_n1 = i_n + 1;
            
            % Compute interpolation weights
            log_T0 = obj.log_temp(i_T0);
            log_T1 = obj.log_temp(i_T1);
            log_n0 = obj.log_dens(i_n0);
            log_n1 = obj.log_dens(i_n1);
            
            dT = log_T1 - log_T0;
            dn = log_n1 - log_n0;
            
            w_T = zeros(size(dT));
            w_n = zeros(size(dn));
            w_T(dT ~= 0) = (log_T(dT ~= 0) - log_T0(dT ~= 0)) ./ dT(dT ~= 0);
            w_n(dn ~= 0) = (log_n(dn ~= 0) - log_n0(dn ~= 0)) ./ dn(dn ~= 0);
            
            % Bilinear interpolation weights
            w00 = (1 - w_T) .* (1 - w_n);
            w01 = (1 - w_T) .* w_n;
            w10 = w_T .* (1 - w_n);
            w11 = w_T .* w_n;
            
            % Extract corner values and interpolate
            E_00 = zeros(n_lines, n_valid);
            E_01 = zeros(n_lines, n_valid);
            E_10 = zeros(n_lines, n_valid);
            E_11 = zeros(n_lines, n_valid);
            
            for j = 1:n_valid
                E_00(:, j) = emiss_table(:, i_T0(j), i_n0(j));
                E_01(:, j) = emiss_table(:, i_T0(j), i_n1(j));
                E_10(:, j) = emiss_table(:, i_T1(j), i_n0(j));
                E_11(:, j) = emiss_table(:, i_T1(j), i_n1(j));
            end
            
            % Apply bilinear interpolation
            emiss_valid = w00' .* E_00 + w01' .* E_01 + w10' .* E_10 + w11' .* E_11;
            emission_rates(in_bounds, :) = emiss_valid';
        end
        
        function [wavelengths, emission_rates] = interpolate_emission_rates_double_vectorized(obj, ...
                Tec_arr, Teh_arr, ne_arr, feh_arr, species_key, min_wav, max_wav)
            % Interpolate per-ion photon emission rates for double Maxwellian.
            %
            % Uses vectorized quadrilinear interpolation in log space.
            % Falls back to single Maxwellian for points where feh or Teh is too small.
            %
            % Parameters
            % ----------
            % Tec_arr : array, shape (n_los,)
            %     Core electron temperatures [eV]
            % Teh_arr : array, shape (n_los,)
            %     Hot electron temperatures [eV]
            % ne_arr : array, shape (n_los,)
            %     Total electron densities [cm^-3]
            % feh_arr : array, shape (n_los,)
            %     Hot electron fractions
            % species_key : char
            %     Species identifier
            % min_wav : float, optional
            %     Minimum wavelength [Angstroms]
            % max_wav : float, optional
            %     Maximum wavelength [Angstroms]
            %
            % Returns
            % -------
            % wavelengths : array, shape (n_lines,)
            %     Emission line wavelengths [Angstroms]
            % emission_rates : array, shape (n_los, n_lines)
            %     Per-ion photon emission rates [photons s^-1 ion^-1]
            
            if nargin < 7, min_wav = 550.0; end
            if nargin < 8, max_wav = 2100.0; end
            
            Tec_arr = Tec_arr(:);
            Teh_arr = Teh_arr(:);
            ne_arr = ne_arr(:);
            feh_arr = feh_arr(:);
            
            % Fall back to single Maxwellian if double not loaded
            if ~obj.double_maxwellian_loaded
                [wavelengths, emission_rates] = obj.interpolate_emission_rates_single_vectorized(...
                    Tec_arr, ne_arr, species_key, min_wav, max_wav);
                return;
            end
            
            if ~isKey(obj.wavelengths_double, species_key)
                wavelengths = [];
                emission_rates = zeros(length(Tec_arr), 0);
                return;
            end
            
            % Get wavelengths and filter
            wav_all = obj.wavelengths_double(species_key);
            wav_mask = (wav_all >= min_wav) & (wav_all <= max_wav);
            wavelengths = wav_all(wav_mask);
            
            if isempty(wavelengths)
                wavelengths = [];
                emission_rates = zeros(length(Tec_arr), 0);
                return;
            end
            
            % Get emission table: shape (n_lines, n_Tc, n_Th, n_n, n_feh)
            emiss_table_all = obj.emissivities_double(species_key);
            emiss_table = emiss_table_all(wav_mask, :, :, :, :);
            
            n_los = length(Tec_arr);
            n_lines = length(wavelengths);
            emission_rates = zeros(n_los, n_lines);
            
            % Determine which points use single vs double Maxwellian
            use_single = (feh_arr < obj.feh_arr(1)) | (Teh_arr < obj.teh_arr(1));
            use_double = ~use_single;
            
            % Handle single Maxwellian fallback
            if any(use_single) && obj.single_maxwellian_loaded
                [~, emiss_s] = obj.interpolate_emission_rates_single_vectorized(...
                    Tec_arr(use_single), ne_arr(use_single), species_key, min_wav, max_wav);
                if ~isempty(emiss_s)
                    emission_rates(use_single, :) = emiss_s;
                end
            end
            
            if ~any(use_double)
                return;
            end
            
            % Check valid double Maxwellian points
            valid = use_double & ...
                   (Tec_arr > 0) & (Teh_arr > 0) & (ne_arr > 0) & ...
                   (feh_arr >= 0) & (feh_arr <= 1) & ...
                   isfinite(Tec_arr) & isfinite(Teh_arr) & isfinite(ne_arr) & isfinite(feh_arr);
            
            in_bounds = valid & ...
                       (Tec_arr >= obj.tec_arr(1)) & (Tec_arr <= obj.tec_arr(end)) & ...
                       (Teh_arr >= obj.teh_arr(1)) & (Teh_arr <= obj.teh_arr(end)) & ...
                       (ne_arr >= obj.ne_arr(1)) & (ne_arr <= obj.ne_arr(end)) & ...
                       (feh_arr >= obj.feh_arr(1)) & (feh_arr <= obj.feh_arr(end));
            
            if ~any(in_bounds)
                return;
            end
            
            % Convert to log space
            log_Tec = log10(Tec_arr(in_bounds));
            log_Teh = log10(Teh_arr(in_bounds));
            log_ne_v = log10(ne_arr(in_bounds));
            log_feh_v = log10(feh_arr(in_bounds));
            n_valid = sum(in_bounds);
            
            % Find bracketing indices
            i_Tec = obj.searchsorted_clip(obj.log_tec, log_Tec, 1, length(obj.log_tec) - 1);
            i_Teh = obj.searchsorted_clip(obj.log_teh, log_Teh, 1, length(obj.log_teh) - 1);
            i_ne = obj.searchsorted_clip(obj.log_ne, log_ne_v, 1, length(obj.log_ne) - 1);
            i_feh = obj.searchsorted_clip(obj.log_feh, log_feh_v, 1, length(obj.log_feh) - 1);
            
            i0_Tec = i_Tec; i1_Tec = i_Tec + 1;
            i0_Teh = i_Teh; i1_Teh = i_Teh + 1;
            i0_ne = i_ne; i1_ne = i_ne + 1;
            i0_feh = i_feh; i1_feh = i_feh + 1;
            
            % Compute interpolation weights
            dTec = obj.log_tec(i1_Tec) - obj.log_tec(i0_Tec);
            dTeh = obj.log_teh(i1_Teh) - obj.log_teh(i0_Teh);
            dne = obj.log_ne(i1_ne) - obj.log_ne(i0_ne);
            dfeh = obj.log_feh(i1_feh) - obj.log_feh(i0_feh);
            
            w_Tec = zeros(size(dTec)); w_Teh = zeros(size(dTeh));
            w_ne = zeros(size(dne)); w_feh = zeros(size(dfeh));
            
            w_Tec(dTec ~= 0) = (log_Tec(dTec ~= 0) - obj.log_tec(i0_Tec(dTec ~= 0))) ./ dTec(dTec ~= 0);
            w_Teh(dTeh ~= 0) = (log_Teh(dTeh ~= 0) - obj.log_teh(i0_Teh(dTeh ~= 0))) ./ dTeh(dTeh ~= 0);
            w_ne(dne ~= 0) = (log_ne_v(dne ~= 0) - obj.log_ne(i0_ne(dne ~= 0))) ./ dne(dne ~= 0);
            w_feh(dfeh ~= 0) = (log_feh_v(dfeh ~= 0) - obj.log_feh(i0_feh(dfeh ~= 0))) ./ dfeh(dfeh ~= 0);
            
            w0_Tec = 1 - w_Tec; w1_Tec = w_Tec;
            w0_Teh = 1 - w_Teh; w1_Teh = w_Teh;
            w0_ne = 1 - w_ne; w1_ne = w_ne;
            w0_feh = 1 - w_feh; w1_feh = w_feh;
            
            % Quadrilinear interpolation over 16 corners of 4D hypercube
            emiss_valid = zeros(n_lines, n_valid);
            
            for b_Tec = 0:1
                for b_Teh = 0:1
                    for b_ne = 0:1
                        for b_feh = 0:1
                            if b_Tec == 0, idx_Tec = i0_Tec; wt_Tec = w0_Tec;
                            else,          idx_Tec = i1_Tec; wt_Tec = w1_Tec; end
                            if b_Teh == 0, idx_Teh = i0_Teh; wt_Teh = w0_Teh;
                            else,          idx_Teh = i1_Teh; wt_Teh = w1_Teh; end
                            if b_ne == 0,  idx_ne = i0_ne;   wt_ne = w0_ne;
                            else,          idx_ne = i1_ne;   wt_ne = w1_ne; end
                            if b_feh == 0, idx_feh = i0_feh; wt_feh = w0_feh;
                            else,          idx_feh = i1_feh; wt_feh = w1_feh; end
                            
                            weight = wt_Tec .* wt_Teh .* wt_ne .* wt_feh;
                            
                            corner_vals = zeros(n_lines, n_valid);
                            for j = 1:n_valid
                                corner_vals(:, j) = emiss_table(:, idx_Tec(j), idx_Teh(j), idx_ne(j), idx_feh(j));
                            end
                            
                            emiss_valid = emiss_valid + weight' .* corner_vals;
                        end
                    end
                end
            end
            
            emission_rates(in_bounds, :) = emiss_valid';
        end
        
        function [wavelengths, brightnesses] = integrate_species_emission_single(obj, ...
                s_values, positions, species_key, Te_los, ne_los, n_ion_los, min_wav, max_wav)
            % Integrate emission for a single species along line of sight (single Maxwellian).
            %
            % Parameters
            % ----------
            % s_values : array
            %     Distances along ray [R_J]
            % positions : array
            %     Positions along ray [R_J]
            % species_key : char
            %     Species identifier
            % Te_los : array
            %     Electron temperature along LOS [eV]
            % ne_los : array
            %     Electron density along LOS [cm^-3]
            % n_ion_los : array
            %     Ion density along LOS [cm^-3]
            % min_wav, max_wav : float
            %     Wavelength range [Angstroms]
            %
            % Returns
            % -------
            % wavelengths : array
            %     Emission line wavelengths [Angstroms]
            % brightnesses : array
            %     Line-integrated brightnesses [Rayleighs]
            
            if nargin < 8, min_wav = 550.0; end
            if nargin < 9, max_wav = 2100.0; end
            
            [wavelengths, emission_rates] = obj.interpolate_emission_rates_single_vectorized(...
                Te_los, ne_los, species_key, min_wav, max_wav);
            
            if isempty(wavelengths)
                wavelengths = [];
                brightnesses = [];
                return;
            end
            
            n_ion_los = n_ion_los(:);
            integrand = emission_rates .* n_ion_los;
            
            n_lines = length(wavelengths);
            brightnesses = zeros(n_lines, 1);
            
            for i = 1:n_lines
                brightnesses(i) = obj.RAYLEIGH_FACTOR * obj.simpson_integrate(integrand(:, i), s_values);
            end
        end
        
        function [wavelengths, brightnesses] = integrate_species_emission_double(obj, ...
                s_values, positions, species_key, Tec_los, Teh_los, ne_los, feh_los, n_ion_los, min_wav, max_wav)
            % Integrate emission for a single species along line of sight (double Maxwellian).
            
            if nargin < 10, min_wav = 550.0; end
            if nargin < 11, max_wav = 2100.0; end
            
            [wavelengths, emission_rates] = obj.interpolate_emission_rates_double_vectorized(...
                Tec_los, Teh_los, ne_los, feh_los, species_key, min_wav, max_wav);
            
            if isempty(wavelengths)
                wavelengths = [];
                brightnesses = [];
                return;
            end
            
            n_ion_los = n_ion_los(:);
            integrand = emission_rates .* n_ion_los;
            
            n_lines = length(wavelengths);
            brightnesses = zeros(n_lines, 1);
            
            for i = 1:n_lines
                brightnesses(i) = obj.RAYLEIGH_FACTOR * obj.simpson_integrate(integrand(:, i), s_values);
            end
        end
        
        function spectrum = convolve_spectrum_erf(~, wavelength_grid, bin_width, ...
                                                   line_wavelengths, line_brightnesses, fwhm)
            % Convolve discrete emission lines with instrument response using ERF.
            %
            % Uses Error Function (ERF) formulation for exact integration of
            % Gaussian line profiles over finite wavelength bins.
            %
            % Parameters
            % ----------
            % wavelength_grid : array
            %     Output wavelength grid (bin centers) [Angstroms]
            % bin_width : float
            %     Width of each wavelength bin [Angstroms]
            % line_wavelengths : array
            %     Discrete emission line wavelengths [Angstroms]
            % line_brightnesses : array
            %     Discrete emission line brightnesses [Rayleighs]
            % fwhm : float, optional
            %     Full width at half maximum of instrument response [Angstroms]
            %
            % Returns
            % -------
            % spectrum : array
            %     Convolved spectrum [Rayleighs/Angstrom]
            
            if nargin < 6, fwhm = 6.0; end
            
            wavelength_grid = wavelength_grid(:);
            
            if isempty(line_wavelengths) || isempty(line_brightnesses)
                spectrum = zeros(size(wavelength_grid));
                return;
            end
            
            % Convert FWHM to Gaussian sigma
            sigma = fwhm / (2 * sqrt(2 * log(2)));
            sigma_sqrt2 = sigma * sqrt(2);
            
            % Calculate bin edges
            bin_lo = wavelength_grid - bin_width / 2;
            bin_hi = wavelength_grid + bin_width / 2;
            
            line_wav = line_wavelengths(:);
            line_bright = line_brightnesses(:);
            
            % Calculate ERF arguments: (n_lines, n_bins)
            erf_arg_lo = (bin_lo' - line_wav) / sigma_sqrt2;
            erf_arg_hi = (bin_hi' - line_wav) / sigma_sqrt2;
            
            % Calculate bin contributions
            bin_contrib = 0.5 * (erf(erf_arg_hi) - erf(erf_arg_lo));
            weighted_contrib = line_bright .* bin_contrib;
            
            % Sum and convert to R/A
            spectrum = sum(weighted_contrib, 1)' / bin_width;
        end
        
        function [wave_bins, spectrum, line_list] = calculate_spectrum_single(obj, ...
                slit_pos_vec, norm_vec, wavelength_range, bin_width, fwhm, ds)
            % Calculate LOS-integrated spectrum for single Maxwellian distribution.
            %
            % Parameters
            % ----------
            % slit_pos_vec : array, shape (3,)
            %     Starting position [x, y, z] in R_J
            % norm_vec : array, shape (3,)
            %     Direction vector (will be normalized)
            % wavelength_range : array, shape (2,), optional
            %     [min_wav, max_wav] in Angstroms (default: [550, 2100])
            % bin_width : float, optional
            %     Spectral bin width [Angstroms] (default: 1.0)
            % fwhm : float, optional
            %     Instrumental FWHM [Angstroms] (default: 6.0)
            % ds : float, optional
            %     Integration step size [R_J] (default: 0.1)
            %
            % Returns
            % -------
            % wave_bins : array
            %     Wavelength bin centers [Angstroms]
            % spectrum : array
            %     Convolved spectrum [Rayleighs/Angstrom]
            % line_list : cell array
            %     List of {wavelength, brightness, species} for each line
            
            if nargin < 4 || isempty(wavelength_range), wavelength_range = [550, 2100]; end
            if nargin < 5 || isempty(bin_width), bin_width = 1.0; end
            if nargin < 6 || isempty(fwhm), fwhm = 6.0; end
            if nargin < 7 || isempty(ds), ds = 0.1; end
            
            fprintf('Calculating single Maxwellian spectrum:\n');
            fprintf('  LOS start: [%.1f, %.1f, %.1f] R_J\n', slit_pos_vec(1), slit_pos_vec(2), slit_pos_vec(3));
            fprintf('  Direction: [%.2f, %.2f, %.2f]\n', norm_vec(1), norm_vec(2), norm_vec(3));
            
            min_wav = wavelength_range(1);
            max_wav = wavelength_range(2);
            
            % Create wavelength grid
            n_wav = floor((max_wav - min_wav) / bin_width) + 1;
            wave_bins = linspace(min_wav, max_wav, n_wav)';
            
            % Trace ray through plasma
            [s_values, positions] = obj.trace_ray(slit_pos_vec, norm_vec, ds);
            
            % Get plasma parameters along LOS
            ne_los = obj.interp_nec(positions(:,1), positions(:,2), positions(:,3));
            Te_los = obj.interp_Tec(positions(:,1), positions(:,2), positions(:,3));
            
            % Clean up NaN values
            ne_los(~isfinite(ne_los) | ne_los <= 0) = 0;
            Te_los(~isfinite(Te_los) | Te_los <= 0) = 0;
            
            % Process each species
            all_wavelengths = [];
            all_brightnesses = [];
            all_species = {};
            
            for k = 1:length(obj.default_species)
                species_key = obj.default_species{k};
                
                if ~isKey(obj.ion_interp_map, species_key)
                    continue;
                end
                
                % Get ion density along LOS
                ion_interp = obj.ion_interp_map(species_key);
                n_ion_los = ion_interp(positions(:,1), positions(:,2), positions(:,3));
                n_ion_los(~isfinite(n_ion_los) | n_ion_los <= 0) = 0;
                
                % Skip if no valid plasma
                if ~any((ne_los > 0) & (Te_los > 0) & (n_ion_los > 0))
                    continue;
                end
                
                % Integrate emission
                [wavelengths, brightnesses] = obj.integrate_species_emission_single(...
                    s_values, positions, species_key, Te_los, ne_los, n_ion_los, min_wav, max_wav);
                
                if ~isempty(wavelengths)
                    all_wavelengths = [all_wavelengths; wavelengths(:)];
                    all_brightnesses = [all_brightnesses; brightnesses(:)];
                    all_species = [all_species; repmat({species_key}, length(wavelengths), 1)];
                end
            end
            
            % Create line list
            line_list = {};
            for i = 1:length(all_wavelengths)
                if all_brightnesses(i) > 1e-10
                    line_list{end+1} = {all_wavelengths(i), all_brightnesses(i), all_species{i}};
                end
            end
            
            fprintf('  Processed %d lines, %d with non-zero brightness\n', ...
                    length(all_wavelengths), length(line_list));
            
            % Convolve with instrument response
            spectrum = obj.convolve_spectrum_erf(wave_bins, bin_width, all_wavelengths, all_brightnesses, fwhm);
        end
        
        function [wave_bins, spectrum, line_list] = calculate_spectrum_double(obj, ...
                slit_pos_vec, norm_vec, wavelength_range, bin_width, fwhm, ds)
            % Calculate LOS-integrated spectrum for double Maxwellian distribution.
            %
            % Same interface as calculate_spectrum_single but uses double Maxwellian
            % emission tables which account for hot electron population effects.
            
            if ~obj.double_maxwellian_loaded
                error('Double Maxwellian tables not loaded.');
            end
            
            if nargin < 4 || isempty(wavelength_range), wavelength_range = [550, 2100]; end
            if nargin < 5 || isempty(bin_width), bin_width = 1.0; end
            if nargin < 6 || isempty(fwhm), fwhm = 6.0; end
            if nargin < 7 || isempty(ds), ds = 0.1; end
            
            fprintf('Calculating double Maxwellian spectrum:\n');
            fprintf('  LOS start: [%.1f, %.1f, %.1f] R_J\n', slit_pos_vec(1), slit_pos_vec(2), slit_pos_vec(3));
            fprintf('  Direction: [%.2f, %.2f, %.2f]\n', norm_vec(1), norm_vec(2), norm_vec(3));
            
            min_wav = wavelength_range(1);
            max_wav = wavelength_range(2);
            
            % Create wavelength grid
            n_wav = floor((max_wav - min_wav) / bin_width) + 1;
            wave_bins = linspace(min_wav, max_wav, n_wav)';
            
            % Trace ray
            [s_values, positions] = obj.trace_ray(slit_pos_vec, norm_vec, ds);
            
            % Get plasma parameters (use total ne for double Maxwellian)
            ne_los = obj.interp_ne_total(positions(:,1), positions(:,2), positions(:,3));
            feh_los = obj.interp_feh(positions(:,1), positions(:,2), positions(:,3));
            Tec_los = obj.interp_Tec(positions(:,1), positions(:,2), positions(:,3));
            Teh_los = obj.interp_Teh(positions(:,1), positions(:,2), positions(:,3));
            
            % Clean up
            ne_los(~isfinite(ne_los) | ne_los <= 0) = 0;
            feh_los(~isfinite(feh_los) | feh_los < 0) = 0;
            Tec_los(~isfinite(Tec_los) | Tec_los <= 0) = 0;
            Teh_los(~isfinite(Teh_los) | Teh_los <= 0) = 0;
            
            % Process each species
            all_wavelengths = [];
            all_brightnesses = [];
            all_species = {};
            
            for k = 1:length(obj.default_species)
                species_key = obj.default_species{k};
                
                if ~isKey(obj.ion_interp_map, species_key)
                    continue;
                end
                
                ion_interp = obj.ion_interp_map(species_key);
                n_ion_los = ion_interp(positions(:,1), positions(:,2), positions(:,3));
                n_ion_los(~isfinite(n_ion_los) | n_ion_los <= 0) = 0;
                
                if ~any((ne_los > 0) & (Tec_los > 0) & (n_ion_los > 0))
                    continue;
                end
                
                [wavelengths, brightnesses] = obj.integrate_species_emission_double(...
                    s_values, positions, species_key, Tec_los, Teh_los, ne_los, feh_los, n_ion_los, min_wav, max_wav);
                
                if ~isempty(wavelengths)
                    all_wavelengths = [all_wavelengths; wavelengths(:)];
                    all_brightnesses = [all_brightnesses; brightnesses(:)];
                    all_species = [all_species; repmat({species_key}, length(wavelengths), 1)];
                end
            end
            
            % Create line list
            line_list = {};
            for i = 1:length(all_wavelengths)
                if all_brightnesses(i) > 1e-10
                    line_list{end+1} = {all_wavelengths(i), all_brightnesses(i), all_species{i}};
                end
            end
            
            fprintf('  Processed %d lines, %d with non-zero brightness\n', ...
                    length(all_wavelengths), length(line_list));
            
            spectrum = obj.convolve_spectrum_erf(wave_bins, bin_width, all_wavelengths, all_brightnesses, fwhm);
        end
        
        function fig = plot_spectrum(obj, wave_bins, spectrum, line_list, varargin)
            % Plot the calculated emission spectrum with optional line markers.
            %
            % Parameters
            % ----------
            % wave_bins : array
            %     Wavelength bin centers [Angstroms]
            % spectrum : array
            %     Spectrum [Rayleighs/Angstrom]
            % line_list : cell array
            %     List of {wavelength, brightness, species} tuples
            % 'title' : char, optional
            %     Plot title
            % 'color' : color specification, optional
            %     Line color (supports 'C0', 'C1' for matplotlib colors)
            %
            % Returns
            % -------
            % fig : figure handle
            
            p = inputParser;
            addParameter(p, 'title', 'Jovian Plasma Torus Emission Spectrum');
            addParameter(p, 'color', [0.1216, 0.4667, 0.7059]);
            parse(p, varargin{:});
            
            plot_title = p.Results.title;
            plot_color = p.Results.color;
            
            % Handle matplotlib color codes
            if ischar(plot_color) || isstring(plot_color)
                if strcmp(plot_color, 'C0')
                    plot_color = [0.1216, 0.4667, 0.7059];  % matplotlib C0 blue
                elseif strcmp(plot_color, 'C1')
                    plot_color = [1.0, 0.4980, 0.0549];     % matplotlib C1 orange
                end
            end
            
            fig = figure('Position', [100, 100, 1200, 600]);
            plot(wave_bins, spectrum, 'Color', plot_color, 'LineWidth', 1.0);
            hold on;
            
            % Mark major emission lines
            if ~isempty(line_list) && length(line_list) > 0
                brightnesses = cellfun(@(x) x{2}, line_list);
                [~, sort_idx] = sort(brightnesses, 'descend');
                
                for i = 1:min(20, length(sort_idx))
                    idx = sort_idx(i);
                    wave = line_list{idx}{1};
                    xline(wave, '--', 'Color', [0.5, 0.5, 0.5], 'Alpha', 0.2, 'LineWidth', 0.5);
                    
                    if i <= 5
                        ion = line_list{idx}{3};
                        if isKey(obj.ion_names, ion)
                            ion_label = obj.ion_names(ion);
                        else
                            ion_label = ion;
                        end
                        ylims = ylim;
                        y_pos = ylims(2) * 0.95;
                        if y_pos <= 0, y_pos = 1.0; end
                        text(wave, y_pos, sprintf('%s\n%.1fA', ion_label, wave), ...
                             'Rotation', 90, 'VerticalAlignment', 'top', ...
                             'HorizontalAlignment', 'right', 'FontSize', 7);
                    end
                end
            end
            
            xlabel('Wavelength [A]', 'FontSize', 12);
            ylabel('Brightness [Rayleighs/A]', 'FontSize', 12);
            title(plot_title, 'FontSize', 14);
            grid on;
            ax = gca;
            ax.GridAlpha = 0.3;
            xlim([wave_bins(1), wave_bins(end)]);
            ylim([0, inf]);
            
            % Add total brightness annotation
            total = obj.simpson_integrate(spectrum, wave_bins);
            annotation('textbox', [0.02, 0.88, 0.15, 0.1], ...
                      'String', sprintf('Total: %.1f R', total), ...
                      'FitBoxToText', 'on', ...
                      'BackgroundColor', 'white', ...
                      'EdgeColor', 'black', ...
                      'FontSize', 10);
            hold off;
        end
        
        function idx = searchsorted_clip(~, arr, values, min_idx, max_idx)
            % Find bracketing indices (equivalent to np.searchsorted with clip).
            %
            % Returns index i such that arr(i) <= value < arr(i+1), clipped to range.
            
            arr = arr(:);
            values = values(:);
            n = length(arr);
            m = length(values);
            idx = zeros(m, 1);
            
            for i = 1:m
                val = values(i);
                lo = 1; hi = n + 1;
                while lo < hi
                    mid = floor((lo + hi) / 2);
                    if arr(mid) < val
                        lo = mid + 1;
                    else
                        hi = mid;
                    end
                end
                % Clip to valid range
                idx(i) = max(min_idx, min(lo - 1, max_idx));
            end
        end
        
        function result = simpson_integrate(~, y, x)
            % Simpson's rule integration (matches scipy.integrate.simpson).
            %
            % Parameters
            % ----------
            % y : array
            %     Function values
            % x : array
            %     Sample points (must be same length as y)
            %
            % Returns
            % -------
            % result : float
            %     Integral approximation
            
            y = y(:);
            x = x(:);
            n = length(y);
            
            if n < 2
                result = 0;
                return;
            elseif n == 2
                result = 0.5 * (y(1) + y(2)) * (x(2) - x(1));
                return;
            elseif n == 3
                h = (x(3) - x(1)) / 2;
                result = h / 3 * (y(1) + 4*y(2) + y(3));
                return;
            end
            
            result = 0;
            if mod(n, 2) == 1
                % Odd number of points - pure Simpson's rule
                for i = 1:2:(n-2)
                    h = (x(i+2) - x(i)) / 2;
                    result = result + h / 3 * (y(i) + 4*y(i+1) + y(i+2));
                end
            else
                % Even number - Simpson's rule + trapezoidal for last interval
                for i = 1:2:(n-3)
                    h = (x(i+2) - x(i)) / 2;
                    result = result + h / 3 * (y(i) + 4*y(i+1) + y(i+2));
                end
                result = result + 0.5 * (y(n-1) + y(n)) * (x(n) - x(n-1));
            end
        end
    end
end
