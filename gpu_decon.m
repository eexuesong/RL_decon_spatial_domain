function output_image = gpu_decon(input_image, FWHM, iterations)
% Updated gpu_decon 4/7/2017
% Features to improve speed and memory utilization
% 1. 3 convolutions with separated Gaussian kernels in X, Y, Z instead of imagaussfilt3.
% This reduces array replication and allows 'in-place' writing - improving memory usage.
% 2. Traffic to/from GPU as uint16 to improve speed

%{
Gaussian function:
f(x) = exp(- x^2 / sigma^2 / 2) / (sigma * sqrt(2 * pi))

The relationship between FWHM and the standard deviation is:
FWHM = 2 * sqrt(2 * ln(2)) * sigma ~= 2.355 * sigma
sigma = FWHM / 2.355
%}

%{
Richardson–Lucy deconvolution:
(1) Blurred_estimate = conv(Estimate(t), PSF)
Note: Estimate(t) means current estimate.

(2) Ratio = Observed_image ./ Blurred_estimate
Note: Observed_image is input_image here.

(3) Correction = conv(Ratio, PSF_flip)

(4) Estimate(t + 1) = Estimate(t) .* Correction
Note: Estimate(t + 1) means output estimate.
%}


% FWHM = [4, 4, 4];   % FWHM pixels values in X, Y, and Z dimensions.
sigma = FWHM / 2.355;   % [1.6985, 1.6985, 1.6985]
% iterations  = 10;

%%%%%%%%%%%%%%%%%% Creating the gaussian kernels first %%%%%%%%%%%%%%%%%%%%
kernal_size = 2 * ceil(2 * sigma) + 1;
% The default filter size is 2*ceil(2*sigma)+1
% About kernal size:
% The Gaussian is infinite in size, but it becomes nearly zero very quickly, and we can truncate it without too much loss.
% The calculation 2 * ceil(2 * sigma) + 1 takes the central portion of the Gaussian of at least four sigma, two sigma to either side.
% The ceiling operation is the “at least”, it needs to be an integer size of course. 
% The +1 accounts for the central pixel.
% This equation always produces an odd size kernel, so it can be symmetric around the origin.
% e.g. kernal_size = [9, 9, 9]

% However, if 2 sigma is quite small for a Gaussian filter 
% (it cuts off too much of the bell shape, affecting some of the good qualities of the filter),
% you can try 3 sigma to either side: 2 * ceil(3 * sigma) + 1.

filter_radius = (kernal_size - 1) / 2;  % [4, 4, 4] 

X = (-filter_radius(1) : filter_radius(1))';    % column-wise [-4; -3; -2; -1; 0; 1; 2; 3; 4]
Y = (-filter_radius(2) : filter_radius(2))';
Z = (-filter_radius(3) : filter_radius(3))';
argX = (X .* X) / (sigma(1) * sigma(1));        % "x^2 / sigma^2" part in Gaussian function: [5.5460; 3.1196; 1.3865; 0.3466; 0; 0.3466; 1.3865; 3.1196; 5.5460]
argY = (Y .* Y) / (sigma(2) * sigma(2));
argZ = (Z .* Z) / (sigma(3) * sigma(3));
psf_X = exp(-argX / 2);                         % "exp(- x^2 / sigma^2 / 2)" part in Gaussian function: [0.0625; 0.2102; 0.4999; 0.8409; 1.0000; 0.8409; 0.4999; 0.2102; 0.0625]
psf_Y = exp(-argY / 2);
psf_Z = exp(-argZ / 2);

%%%%%%%%%%%%%%%%%%%%% Suppress near-zero components %%%%%%%%%%%%%%%%%%%%%%%
% In practice, has no effect. Just keep.
psf_X(psf_X < eps * max(psf_X(:))) = 0;
psf_Y(psf_Y < eps * max(psf_Y(:))) = 0;
psf_Z(psf_Z < eps * max(psf_Z(:))) = 0;

%%%%%%%%%%%%%%%%%%%%%% PSF normalization and reshape %%%%%%%%%%%%%%%%%%%%%%
% X dimension
sum_psf_X = sum(psf_X(:));
psf_X = single(psf_X ./ sum_psf_X);
psf_X = reshape(psf_X, [1, kernal_size(1)]);    % [0.014779855, 0.049722526, 0.11827642, 0.19893225, 0.23657791, 0.19893225, 0.11827642, 0.049722526, 0.014779855]

% Y dimension
sumg_psf_Y = sum(psf_Y(:));
psf_Y = single(psf_Y ./ sumg_psf_Y);            % [0.014779855; 0.049722526; 0.11827642; 0.19893225; 0.23657791; 0.19893225; 0.11827642; 0.049722526; 0.014779855]

% Z dimension
sum_psf_Z = sum(psf_Z(:));
psf_Z = single(psf_Z ./ sum_psf_Z);
psf_Z = reshape(psf_Z, [1, 1, kernal_size(3)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pad input image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ny, nx, nz] = size(input_image);       % e.g. [512, 512, 100]
npad = ceil(3 * sigma);                 % [6, 6, 6]
gpu_Estimate = gpuArray(single(input_image));
if nz == 1
    % 2D slice
    gpu_Estimate = padarray(gpu_Estimate, [npad(2), npad(1)], 'replicate', 'both');
else
    % 3D stack
    gpu_Estimate = padarray(gpu_Estimate, [npad(2), npad(1), npad(3)], 'replicate', 'both');
end


%%%%%%%%%%%%%%%%%%%%%%%%%% Normalization factor %%%%%%%%%%%%%%%%%%%%%%%%%%%
Ratio = ones(size(gpu_Estimate), 'single', 'gpuArray');
gpu_Normalization_factor = convn(Ratio, psf_X, 'same');
gpu_Normalization_factor = convn(gpu_Normalization_factor, psf_Y, 'same');
if nz > 1
    % 3D stack
    gpu_Normalization_factor = convn(gpu_Normalization_factor, psf_Z, 'same');
end
gpu_Normalization_factor(gpu_Normalization_factor <= 0) = eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Back projector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gpu_input_image = gpu_Estimate;

for k = 1: iterations
    Blur = convn(gpu_Estimate, psf_X, 'same');
    Blur = convn(Blur, psf_Y, 'same');
    if nz > 1
        Blur = convn(Blur, psf_Z, 'same');
    end
    % 'same' — Return the central part of the convolution, which is the same size as A in convn(A,B).

    Ratio = gpu_input_image ./ Blur;
    Correction = convn(Ratio, psf_X, 'same');
    Correction = convn(Correction, psf_Y, 'same');
    if nz > 1
        Correction = convn(Correction, psf_Z, 'same');
    end
    % Because psf_X, psf_Y and psf_Z are all symmetric around the origin.
    % Their flip is same as themselves.

    gpu_Estimate = gpu_Estimate .* Correction ./ gpu_Normalization_factor;
end

if nz == 1
    % 2D slice
    gpu_Estimate = gpu_Estimate(npad(2)+1:npad(2)+ny, npad(1)+1:npad(1)+nx);
else
    % 3D stack
    gpu_Estimate = gpu_Estimate(npad(2)+1:npad(2)+ny, npad(1)+1:npad(1)+nx, npad(3)+1:npad(3)+nz);
end
output_image = gather(gpu_Estimate);
% output_image = uint16(65535 * output_image ./ max(output_image, [], 'all'));


