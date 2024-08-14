clear all;

%% X-Y-Z 3D stack
% FWHM for x dimension: 4 pixels (184 nm)
% FHWM for y dimension: 4 pixels (184 nm)
% FHWM for z dimension: 4 pixels (400 nm)
% Background: 100 counts
% Iteration number: 10 times
% decon_gpu_3d(4, 4, 4, 100, 10, "D:\Code\Matlab_Code\RL_deconvolution_spatial_domain\C2-Co-Local__1_MMStack_Default.ome.tif");

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



% FWHM for x-y dimension:
% #1.5 coverslip only
%       1st channel (green): 3.5 pixels (402.5 nm)
%       2nd channel (red):   3.5 pixels (402.5 nm)
% FWHM for z dimension:
%       1st channel (green): 9 pixels (900 nm)
%       2nd channel (red):   9 pixels (900 nm)

% coverslide + #1.5 coverslip
%       1st channel (green): 3 pixels (345 nm)
%       2nd channel (red):   3 pixels (345 nm)
% FWHM for z dimension:
%       1st channel (green): 9.5 pixels (950 nm)
%       2nd channel (red):   9.5 pixels (950 nm)

% Background: 200 counts
% Iteration number: 10 / 30 times
decon_3d([3, 3], [3, 3], [9.5, 9.5], [200, 200], [10, 10], "Z:\shrofflab\Alyssa\080924_Col-76_Dpy-7_Cross\Coverslip + Slide\080924_C+S_Pos_3\080924_C+S_Pos_3_MMStack_Default.ome.tif");