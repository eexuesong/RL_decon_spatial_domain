clear all;

%% X-Y-Z 3D stack
% FWHM for x dimension: 4 pixels (184 nm)
% FHWM for y dimension: 4 pixels (184 nm)
% FHWM for z dimension: 4 pixels (400 nm)
% Iteration number: 10 times
gpu_decon_3d(4, 4, 4, 10, "C2-Co-Local__1_MMStack_Default.ome.tif");

%% X-Y-C-Z 4D stack
% FWHM for x dimension:
%       1st channel (red): 4.5 pixels (207 nm)
%       2nd channel (green): 4 pixels (184 nm)
% FWHM for y dimension:
%       1st channel (red): 4.5 pixels (207 nm)
%       2nd channel (green): 4 pixels (184 nm)
% FWHM for z dimension:
%       1st channel (red): 4.5 pixels (450 nm)
%       2nd channel (green): 4 pixels (400 nm)
% Iteration number:
%       1st channel (red): 10 times
%       2nd channel (gree): 10 times
gpu_decon_3d([4.5, 4], [4.5, 4], [4.5, 4], [10, 10], "Co-Local__1_MMStack_Default.ome.tif");