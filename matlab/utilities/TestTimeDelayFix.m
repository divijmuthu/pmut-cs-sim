% Test Time Delay Fix
% Quick test to verify the time delay calculation is now correct
clearvars;
clc;

%% Test the corrected time delay calculation
fprintf('--- Testing Time Delay Fix ---\n');

% Test parameters
params.c = 343;  % Speed of sound [m/s]
distance_mm = 267;  % Example distance in mm

% OLD (wrong) calculation
time_delay_old = distance_mm / (params.c * 1e-3);
fprintf('OLD calculation: %.6f seconds (WRONG!)\n', time_delay_old);

% NEW (correct) calculation
time_delay_new = distance_mm / (params.c * 1e3);
fprintf('NEW calculation: %.6f seconds (CORRECT!)\n', time_delay_new);

% Convert to microseconds for comparison
time_delay_us = time_delay_new * 1e6;
fprintf('Time delay: %.2f microseconds (reasonable!)\n', time_delay_us);

% Test response calculation
freq = 57000;  % Hz
wavelength = params.c / freq;  % m
wavelength_mm = wavelength * 1000;  % mm

spatial_response = exp(-distance_mm / (wavelength_mm * 2));
phase_factor = cos(2 * pi * distance_mm / wavelength_mm);
response = spatial_response * phase_factor * exp(-time_delay_new * 1e3) * 1e-1;

fprintf('Response calculation:\n');
fprintf('  Distance: %.1f mm\n', distance_mm);
fprintf('  Wavelength: %.1f mm\n', wavelength_mm);
fprintf('  Time delay: %.6f s (%.2f us)\n', time_delay_new, time_delay_us);
fprintf('  Spatial response: %.6f\n', spatial_response);
fprintf('  Phase factor: %.6f\n', phase_factor);
fprintf('  Final response: %.6f\n', response);

if response > 1e-10
    fprintf('✅ Response is now non-zero!\n');
else
    fprintf('❌ Response is still too small\n');
end

fprintf('\n--- Test Complete ---\n'); 