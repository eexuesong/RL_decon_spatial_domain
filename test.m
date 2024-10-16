clear all;

%% X-Y-T 3D stack
% FWHM for x dimension: 4 pixels (184 nm)
% FHWM for y dimension: 4 pixels (184 nm)
% Background: 100 counts
% Iteration number: 10 times
decon_2d("cpu", 4, 4, 100, 10, "D:\Code\Matlab_Code\RL_deconvolution_spatial_domain\ch3_flow-0mlmin_3\ch3_flow-0mlmin_3_MMStack_Default.ome-1.tif");

%% X-Y-Z 3D stack
% FWHM for x dimension: 4 pixels (184 nm)
% FHWM for y dimension: 4 pixels (184 nm)
% FHWM for z dimension: 4 pixels (400 nm)
% Iteration number: 10 times
% decon_3d(4, 4, 4, 100, 10, "D:\Code\Matlab_Code\RL_deconvolution_spatial_domain\C2-Co-Local__1_MMStack_Default.ome.tif");

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
% decon_3d([4.5, 4], [4.5, 4], [4.5, 4], [100, 100], [10, 10], "Co-Local__1_MMStack_Default.ome.tif");