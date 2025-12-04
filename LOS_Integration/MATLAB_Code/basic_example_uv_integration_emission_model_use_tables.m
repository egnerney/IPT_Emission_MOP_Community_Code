function results = basic_example_uv_integration_emission_model_use_tables()
% basic_example_uv_integration_emission_model_use_tables
%
% Example UV Emission Calculation Using Line-of-Sight Integration
% ================================================================
%
% This script demonstrates how to use the Jovian UV emission raytracer to calculate
% synthetic UV spectra through the Io Plasma Torus using proper line-of-sight
% integration with species-by-species interpolation from CHIANTI emission tables.
%
% The example calculates UV emission spectra for test cases with both single and
% double Maxwellian electron distributions:
% 1. Single Maxwellian: equatorial and off-equator lines of sight
% 2. Double Maxwellian: same geometries with hot electron population
%
% All lines of sight start at x=6 R_J and look through the plasma torus in the
% +y direction, sampling the peak emission region near Io's orbit.
%
% PHYSICAL MODEL:
% - Ray tracing through 3D plasma model with trilinear interpolation
% - Vectorized bilinear interpolation of single Maxwellian emission rates (2D tables)
% - Vectorized quadrilinear interpolation of double Maxwellian emission rates (4D tables)
% - Per-ion photon emission rates multiplied by ion density and 1e-6 Rayleigh factor
% - Simpson's rule integration over line of sight
% - Analytic ERF-based Gaussian convolution for instrument response
%
% DIRECTORY STRUCTURE EXPECTED:
% - IPT_Emission_MOP_Community_Code/
%   - LOS_Integration/MATLAB_Code/           (this script and module)
%   - Emiss_tables/                          (CHIANTI emission tables)
%   - 3D_Torus_Model/                        (3D plasma model in HDF5 format)
%
% OUTPUT FILES:
% - spectrum_single_equatorial.png: Single Maxwellian equatorial spectrum
% - spectrum_single_off_equator.png: Single Maxwellian off-equator spectrum
% - spectrum_single_comparison.png: Comparison of equatorial vs off-equator
% - spectrum_double_equatorial.png: Double Maxwellian equatorial spectrum
% - spectrum_double_off_equator.png: Double Maxwellian off-equator spectrum
% - spectrum_single_vs_double.png: Comparison of single vs double Maxwellian
%
% AUTHOR: Edward (Eddie) G. Nerney
% INSTITUTION: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
% DATE: November 2025
% LICENSE: Open source for academic and research use

fprintf('======================================================================\n');
fprintf('Jovian UV Emission Raytracer - Line-of-Sight Integration\n');
fprintf('======================================================================\n');
fprintf('\n');
fprintf('This example demonstrates UV emission calculations through\n');
fprintf('the Io Plasma Torus using proper line-of-sight integration\n');
fprintf('with species-by-species interpolation from CHIANTI tables.\n');
fprintf('\n');

start_time = tic;

% ========================================================================
% INITIALIZE RAYTRACER
% ========================================================================

fprintf('Loading plasma model and emission tables...\n');
try
    raytracer = IPT_emiss_MOP_community_code();
catch ME
    fprintf('\nError: %s\n', ME.message);
    fprintf('\nPlease ensure the following directory structure exists:\n');
    fprintf('  IPT_Emission_MOP_Community_Code/\n');
    fprintf('    - LOS_Integration/MATLAB_Code/ (containing this script)\n');
    fprintf('    - Emiss_tables/CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5\n');
    fprintf('    - Emiss_tables/CHIANTI_11.0.2_emiss_tables_double_maxwellian_16x8x12x10.h5\n');
    fprintf('    - 3D_Torus_Model/jovian_plasma_interpolated_381x381x231.h5\n');
    results = [];
    return;
end

% ========================================================================
% DEFINE OBSERVATION GEOMETRY
% ========================================================================

% Common parameters
wavelength_range = [550, 2100];  % UV range [Angstroms]
bin_width = 1.0;                 % Spectral bin width [Angstroms]
fwhm = 6.0;                      % Instrumental FWHM [Angstroms] - Europa-UVS/JUICE-UVS
ds = 0.1;                        % Integration step size [R_J]

% Test case 1: Equatorial line of sight
slit_pos_vec1 = [6.0, -20.0, 0.0];  % Start at x=6 R_J, equator
norm_vec1 = [0.0, 1.0, 0.0];         % Look in +y direction

% Test case 2: Off-equator line of sight
slit_pos_vec2 = [6.0, -20.0, 0.5];   % 0.5 R_J above equator
norm_vec2 = [0.0, 1.0, 0.0];          % Look in +y direction

% Matplotlib color codes for consistency with Python plots
color_C0 = [0.1216, 0.4667, 0.7059];  % #1f77b4 blue
color_C1 = [1.0, 0.4980, 0.0549];     % #ff7f0e orange

% ========================================================================
% SINGLE MAXWELLIAN CALCULATIONS
% ========================================================================

fprintf('\n======================================================================\n');
fprintf('Test Case 1: Equatorial Line of Sight (Single Maxwellian)\n');
fprintf('----------------------------------------------------------------------\n');
fprintf('Starting position: (%.1f, %.1f, %.1f) R_J\n', slit_pos_vec1(1), slit_pos_vec1(2), slit_pos_vec1(3));
fprintf('Direction vector: (%.1f, %.1f, %.1f)\n', norm_vec1(1), norm_vec1(2), norm_vec1(3));
fprintf('\n');

[wave_bins1_single, spectrum1_single, lines1_single] = raytracer.calculate_spectrum_single(...
    slit_pos_vec1, norm_vec1, wavelength_range, bin_width, fwhm, ds);

% Calculate summary statistics
total_brightness1_single = simpson_integrate_standalone(spectrum1_single, wave_bins1_single);
[~, peak_idx1_single] = max(spectrum1_single);

fprintf('\nSingle Maxwellian Results:\n');
fprintf('  Total brightness: %.1f Rayleighs\n', total_brightness1_single);
fprintf('  Peak brightness: %.2f R/A at %.1f A\n', max(spectrum1_single), wave_bins1_single(peak_idx1_single));

% Test case 2: Off-equator (single Maxwellian)
fprintf('\n======================================================================\n');
fprintf('Test Case 2: Off-Equator Line of Sight (Single Maxwellian)\n');
fprintf('----------------------------------------------------------------------\n');
fprintf('Starting position: (%.1f, %.1f, %.1f) R_J\n', slit_pos_vec2(1), slit_pos_vec2(2), slit_pos_vec2(3));
fprintf('Direction vector: (%.1f, %.1f, %.1f)\n', norm_vec2(1), norm_vec2(2), norm_vec2(3));
fprintf('\n');

[wave_bins2_single, spectrum2_single, lines2_single] = raytracer.calculate_spectrum_single(...
    slit_pos_vec2, norm_vec2, wavelength_range, bin_width, fwhm, ds);

total_brightness2_single = simpson_integrate_standalone(spectrum2_single, wave_bins2_single);
[~, peak_idx2_single] = max(spectrum2_single);

fprintf('\nSingle Maxwellian Results:\n');
fprintf('  Total brightness: %.1f Rayleighs\n', total_brightness2_single);
fprintf('  Peak brightness: %.2f R/A at %.1f A\n', max(spectrum2_single), wave_bins2_single(peak_idx2_single));

% ========================================================================
% DOUBLE MAXWELLIAN CALCULATIONS
% ========================================================================

% Initialize variables
wave_bins1_double = [];
spectrum1_double = [];
lines1_double = {};
wave_bins2_double = [];
spectrum2_double = [];
lines2_double = {};
total_brightness1_double = 0;
total_brightness2_double = 0;
enhancement1 = 1.0;
enhancement2 = 1.0;

if raytracer.double_maxwellian_loaded
    % Test case 3: Equatorial (double Maxwellian)
    fprintf('\n======================================================================\n');
    fprintf('Test Case 3: Equatorial LOS (Double Maxwellian)\n');
    fprintf('----------------------------------------------------------------------\n');
    
    [wave_bins1_double, spectrum1_double, lines1_double] = raytracer.calculate_spectrum_double(...
        slit_pos_vec1, norm_vec1, wavelength_range, bin_width, fwhm, ds);
    
    total_brightness1_double = simpson_integrate_standalone(spectrum1_double, wave_bins1_double);
    [~, peak_idx1_double] = max(spectrum1_double);
    enhancement1 = total_brightness1_double / total_brightness1_single;
    
    fprintf('\nDouble Maxwellian Results:\n');
    fprintf('  Total brightness: %.1f Rayleighs\n', total_brightness1_double);
    fprintf('  Peak brightness: %.2f R/A at %.1f A\n', max(spectrum1_double), wave_bins1_double(peak_idx1_double));
    fprintf('  Enhancement factor: %.3f\n', enhancement1);
    
    % Test case 4: Off-equator (double Maxwellian)
    fprintf('\n======================================================================\n');
    fprintf('Test Case 4: Off-Equator LOS (Double Maxwellian)\n');
    fprintf('----------------------------------------------------------------------\n');
    
    [wave_bins2_double, spectrum2_double, lines2_double] = raytracer.calculate_spectrum_double(...
        slit_pos_vec2, norm_vec2, wavelength_range, bin_width, fwhm, ds);
    
    total_brightness2_double = simpson_integrate_standalone(spectrum2_double, wave_bins2_double);
    [~, peak_idx2_double] = max(spectrum2_double);
    enhancement2 = total_brightness2_double / total_brightness2_single;
    
    fprintf('\nDouble Maxwellian Results:\n');
    fprintf('  Total brightness: %.1f Rayleighs\n', total_brightness2_double);
    fprintf('  Peak brightness: %.2f R/A at %.1f A\n', max(spectrum2_double), wave_bins2_double(peak_idx2_double));
    fprintf('  Enhancement factor: %.3f\n', enhancement2);
end

% ========================================================================
% GENERATE PLOTS
% ========================================================================

fprintf('\n======================================================================\n');
fprintf('Generating Plots\n');
fprintf('----------------------------------------------------------------------\n');

% Plot 1: Single Maxwellian equatorial
fig1 = raytracer.plot_spectrum(wave_bins1_single, spectrum1_single, lines1_single, ...
    'title', 'Single Maxwellian: Equatorial LOS (x=6 R_J, z=0)', 'color', color_C0);
saveas(fig1, 'spectrum_single_equatorial.png');
fprintf('Saved: spectrum_single_equatorial.png\n');

% Plot 2: Single Maxwellian off-equator
fig2 = raytracer.plot_spectrum(wave_bins2_single, spectrum2_single, lines2_single, ...
    'title', 'Single Maxwellian: Off-Equator LOS (x=6 R_J, z=0.5 R_J)', 'color', color_C0);
saveas(fig2, 'spectrum_single_off_equator.png');
fprintf('Saved: spectrum_single_off_equator.png\n');

% Plot 3: Comparison - equatorial vs off-equator (single Maxwellian)
fig3 = figure('Position', [100, 100, 1200, 600]);
plot(wave_bins1_single, spectrum1_single, 'Color', color_C0, 'LineWidth', 1.0, ...
     'DisplayName', 'Equatorial (z=0)');
hold on;
plot(wave_bins2_single, spectrum2_single, 'Color', color_C1, 'LineWidth', 1.0, ...
     'DisplayName', 'Off-equator (z=0.5 R_J)');
hold off;
xlabel('Wavelength [A]', 'FontSize', 12);
ylabel('Brightness [Rayleighs/A]', 'FontSize', 12);
title('Single Maxwellian: Equatorial vs Off-Equator UV Spectra', 'FontSize', 14);
grid on;
ax = gca;
ax.GridAlpha = 0.3;
legend('Location', 'northeast');
xlim(wavelength_range);
ylim([0, inf]);
annotation('textbox', [0.02, 0.88, 0.15, 0.1], ...
          'String', sprintf('Equatorial: %.1f R\nOff-equator: %.1f R', ...
                           total_brightness1_single, total_brightness2_single), ...
          'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 10);
saveas(fig3, 'spectrum_single_comparison.png');
fprintf('Saved: spectrum_single_comparison.png\n');

% Double Maxwellian plots (if available)
if raytracer.double_maxwellian_loaded && ~isempty(spectrum1_double)
    % Plot 4: Double Maxwellian equatorial
    fig4 = raytracer.plot_spectrum(wave_bins1_double, spectrum1_double, lines1_double, ...
        'title', 'Double Maxwellian: Equatorial LOS', 'color', color_C1);
    saveas(fig4, 'spectrum_double_equatorial.png');
    fprintf('Saved: spectrum_double_equatorial.png\n');
    
    % Plot 5: Double Maxwellian off-equator
    fig5 = raytracer.plot_spectrum(wave_bins2_double, spectrum2_double, lines2_double, ...
        'title', 'Double Maxwellian: Off-Equator LOS', 'color', color_C1);
    saveas(fig5, 'spectrum_double_off_equator.png');
    fprintf('Saved: spectrum_double_off_equator.png\n');
    
    % Plot 6: Comparison - single vs double Maxwellian (equatorial)
    fig6 = figure('Position', [100, 100, 1200, 600]);
    plot(wave_bins1_single, spectrum1_single, 'Color', color_C0, 'LineWidth', 1.0, ...
         'DisplayName', 'Single Maxwellian');
    hold on;
    plot(wave_bins1_double, spectrum1_double, 'Color', color_C1, 'LineWidth', 1.0, ...
         'DisplayName', 'Double Maxwellian');
    hold off;
    xlabel('Wavelength [A]', 'FontSize', 12);
    ylabel('Brightness [Rayleighs/A]', 'FontSize', 12);
    title('Equatorial LOS: Single vs Double Maxwellian', 'FontSize', 14);
    grid on;
    ax = gca;
    ax.GridAlpha = 0.3;
    legend('Location', 'northeast');
    xlim(wavelength_range);
    ylim([0, inf]);
    annotation('textbox', [0.02, 0.88, 0.18, 0.1], ...
              'String', sprintf('Single: %.1f R\nDouble: %.1f R\nEnhancement: %.3fx', ...
                               total_brightness1_single, total_brightness1_double, enhancement1), ...
              'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 10);
    saveas(fig6, 'spectrum_single_vs_double.png');
    fprintf('Saved: spectrum_single_vs_double.png\n');
end

% ========================================================================
% STRONGEST EMISSION LINES
% ========================================================================

fprintf('\n======================================================================\n');
fprintf('Strongest UV Emission Lines (Equatorial, Single Maxwellian)\n');
fprintf('----------------------------------------------------------------------\n');

if ~isempty(lines1_single)
    % Sort by brightness
    brightnesses = cellfun(@(x) x{2}, lines1_single);
    [~, sort_idx] = sort(brightnesses, 'descend');
    
    ion_display = containers.Map(...
        {'SP', 'S2P', 'S3P', 'S4P', 'OP', 'O2P'}, ...
        {'S II', 'S III', 'S IV', 'S V', 'O II', 'O III'});
    
    fprintf('%-8s %-12s %-12s\n', 'Ion', 'Wavelength', 'Brightness');
    fprintf('%s\n', repmat('-', 1, 35));
    
    for i = 1:min(15, length(sort_idx))
        idx = sort_idx(i);
        wave = lines1_single{idx}{1};
        brightness = lines1_single{idx}{2};
        ion = lines1_single{idx}{3};
        
        if isKey(ion_display, ion)
            ion_label = ion_display(ion);
        else
            ion_label = ion;
        end
        
        fprintf('%-8s %8.2f A   %8.2f R\n', ion_label, wave, brightness);
    end
end

% ========================================================================
% DIAGNOSTIC INFORMATION
% ========================================================================

fprintf('\n======================================================================\n');
fprintf('Diagnostic Information: Plasma Parameters Along Ray (Equatorial)\n');
fprintf('----------------------------------------------------------------------\n');

% Trace ray and sample plasma parameters
[s_test, pos_test] = raytracer.trace_ray(slit_pos_vec1, norm_vec1, 1.0);

if length(s_test) > 0
    fprintf('\nSampling plasma parameters at key positions:\n');
    fprintf('%-30s %-10s %-12s %-10s\n', 'Position (R_J)', 'rho (R_J)', 'ne (cm^-3)', 'Te (eV)');
    fprintf('%s\n', repmat('-', 1, 65));
    
    % Sample at key points along ray
    if length(s_test) > 20
        sample_indices = [11, floor(length(s_test)/4), floor(length(s_test)/2), ...
                         floor(3*length(s_test)/4), length(s_test)-9];
    else
        sample_indices = [1, floor(length(s_test)/2)+1, length(s_test)];
    end
    
    for idx = sample_indices
        if idx <= size(pos_test, 1)
            pos = pos_test(idx, :);
            ne = raytracer.interp_nec(pos(1), pos(2), pos(3));
            Te = raytracer.interp_Tec(pos(1), pos(2), pos(3));
            if isnan(ne), ne = 0; end
            if isnan(Te), Te = 0; end
            rho = sqrt(pos(1)^2 + pos(2)^2);
            
            pos_str = sprintf('[%6.1f, %6.1f, %6.1f]', pos(1), pos(2), pos(3));
            fprintf('%-30s %8.2f   %10.2f   %8.2f\n', pos_str, rho, ne, Te);
        end
    end
end

% ========================================================================
% PHYSICAL INTERPRETATION
% ========================================================================

fprintf('\n======================================================================\n');
fprintf('Physical Interpretation\n');
fprintf('----------------------------------------------------------------------\n');
fprintf('The plasma torus shows:\n');
fprintf('- Peak emission near 6 R_J (Io''s orbital radius)\n');
fprintf('- Scale height of ~0.5-1 R_J\n');
fprintf('- Dominant emission from S III and O II ions\n');
fprintf('- Temperature ~5-10 eV in the cold torus\n');
fprintf('- Electron density ~100-2000 cm^-3 at peak\n');
if raytracer.double_maxwellian_loaded
    fprintf('- Hot electrons enhance high-excitation lines\n');
    fprintf('- Overall brightness enhancement: %.1f%%\n', (enhancement1-1)*100);
end

% ========================================================================
% TIMING AND SUMMARY
% ========================================================================

elapsed_time = toc(start_time);

fprintf('\n======================================================================\n');
fprintf('UV Emission Calculation Complete\n');
fprintf('======================================================================\n');
fprintf('\nExecution time: %.2f seconds\n', elapsed_time);
fprintf('\nOutput files generated:\n');
fprintf('  - spectrum_single_equatorial.png\n');
fprintf('  - spectrum_single_off_equator.png\n');
fprintf('  - spectrum_single_comparison.png\n');
if raytracer.double_maxwellian_loaded
    fprintf('  - spectrum_double_equatorial.png\n');
    fprintf('  - spectrum_double_off_equator.png\n');
    fprintf('  - spectrum_single_vs_double.png\n');
end

fprintf('\nCalculation methodology:\n');
fprintf('  - Ray tracing through 3D plasma model with trilinear interpolation\n');
fprintf('  - Vectorized bilinear interpolation of single Maxwellian emission rates (2D)\n');
if raytracer.double_maxwellian_loaded
    fprintf('  - Vectorized quadrilinear interpolation of double Maxwellian emission rates (4D)\n');
end
fprintf('  - Per-ion photon emission rates x ion density x 1e-6 -> Rayleighs\n');
fprintf('  - Simpson''s rule for line-of-sight integration\n');
fprintf('  - ERF-based Gaussian convolution for instrument response\n');
fprintf('  - Integration step size: ds = %.2f R_J\n', ds);
fprintf('  - Instrumental FWHM: %.1f A\n', fwhm);

% Return results dictionary for further analysis
results = struct();
results.wave_bins_single = wave_bins1_single;
results.spectrum_single_equatorial = spectrum1_single;
results.spectrum_single_off_equator = spectrum2_single;
results.lines_single_equatorial = lines1_single;
results.lines_single_off_equator = lines2_single;
results.total_brightness_single_equatorial = total_brightness1_single;
results.total_brightness_single_off_equator = total_brightness2_single;
results.raytracer = raytracer;

if raytracer.double_maxwellian_loaded && ~isempty(spectrum1_double)
    results.wave_bins_double = wave_bins1_double;
    results.spectrum_double_equatorial = spectrum1_double;
    results.spectrum_double_off_equator = spectrum2_double;
    results.lines_double_equatorial = lines1_double;
    results.lines_double_off_equator = lines2_double;
    results.total_brightness_double_equatorial = total_brightness1_double;
    results.total_brightness_double_off_equator = total_brightness2_double;
    results.enhancement_equatorial = enhancement1;
    results.enhancement_off_equator = enhancement2;
end

end


function result = simpson_integrate_standalone(y, x)
% Simpson's rule integration (standalone version matching scipy.integrate.simpson).
%
% Parameters
% ----------
% y : array
%     Function values to integrate
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
    % Trapezoidal rule for 2 points
    result = 0.5 * (y(1) + y(2)) * (x(2) - x(1));
    return;
elseif n == 3
    % Pure Simpson's rule for 3 points
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
