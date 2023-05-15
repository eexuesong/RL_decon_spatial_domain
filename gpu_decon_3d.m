function gpu_decon_3d(xres, yres, zres, iteration, filename)
% Reads a image file, processes it, saves an output tif as 'filename_decon.tif'.

% 09/06/17: gpu based serial decon, using gpu_decon. (USE THIS ONE)

% The Filename e should not contain space
% xres yres zres: In pixels, the full width at half maximum of the Gaussian PSF.

% [4, 4, 5]:
% The actual FWHM in X-Y is around 180 nm. Thefore the value = 0.18 / pixel_size = 0.18 / 0.046 = 3.9 pixels.
% X,Y is normally keep constant which is always around 4 pixels.
% The actual FWHM in Z is around 400 nm. Therefore the value = 0.4 /z_step_size = 0.4 / 0.1 = 4 pixels.

% This is the newer version of Decon. Jiji 07-01-2022
% Correct way to use the program for example:
% batch_gpu_decon_3d(4, 4, 4, 10, "D:\Code\Matlab_Code\RL_deconvolution")

% For 60x water 1.2 NA lens FWHM in X-Y is [2.5, 2.5] and in Z equals to 0.5 / z_step_size.
% For ORCA_Quest, pixel size 46 nm. For 100x oil lens The X-Y should be 4 pixels without binning. And 2 with 2x2 binning.

% For 60x water lens, the X-Y should be [3.9, 3.9].
% Z FWHM is measured around 0.58. Therefore, Z value should be 0.58 / z_step_size


%% Read tiff stack file
[image, header] = ImageJ_formatted_TIFF.ReadTifStack(filename);
if isempty(header.channels)
    channels = 1;
else
    channels = header.channels;
    if (channels ~= size(xres, 2)) || (channels ~= size(yres, 2)) || (channels ~= size(zres, 2))
        error("Channel number of image stack is inconsistent with xres, yres, zres or iteration array size.");
    end
end

% 2D slice or 3D stack
if isempty(header.slices)
    slices = 1;
else
    slices = header.slices;
end

output_stack = zeros(size(image), 'uint16');
if (channels ~= 1) && (slices ~= 1)
    % 4D 'cz' stack only
    output_stack = permute(output_stack, [1 2 4 3]); % Swap channels and slices
end

%% RL deconvolution part
g = gpuDevice(1); reset(g);
s = imfinfo(filename);

if s(1).FileSize <= (g.FreeMemory / 4.4)
    disp(['GPU Memory before RL deconvolution: ',num2str(g.FreeMemory / 1024 / 1024 / 1024), ' GB']);
    gpu_flag = 1;
else
    gpu_flag = 0;
end

for channel_index = 1:channels
    if channels == 1
        input_image = image;
    else
        input_image = image(:, :, channel_index, :);
        if slices > 1
            input_image = squeeze(permute(input_image, [1 2 4 3]));    % Swap slices and channels
        end
    end

%     if min(input_image, [], 'all') <= 0
%         input_image = input_image + 0.00001;
%     end

    FWHM = [xres(channel_index), yres(channel_index), zres(channel_index)];
    if gpu_flag
        output_image = gpu_decon(input_image, FWHM, iteration(channel_index));
    else
        output_image = cpu_decon(input_image, FWHM, iteration(channel_index));
    end

    if (header.BitsPerSample == 16)
        if max(output_image, [], 'all') <= 65535
            output_image = uint16(output_image);
        else
            output_image = uint16(65535 * output_image ./ max(output_image, [], 'all'));
        end
    end
    
    if slices == 1
        output_stack(:, :, channel_index) = output_image;
    else
        output_stack(:, :, :, channel_index) = output_image;
    end

    if gpu_flag
        disp(['GPU Memory after RL deconvolution: ',num2str(g.FreeMemory / 1024 / 1024 / 1024), ' GB']);
    end
end


%% Save tiff stack file
[filepath, name, ext] = fileparts(filename);
path_output = strcat(filepath, '\deconvolution_results\');
mkdir(path_output);
filepath_output = strcat(path_output, name, '_decon.tif');

if channels == 1
    if isempty(header.resolution)
        ImageJ_formatted_TIFF.WriteTifStack(output_stack, filepath_output);
    elseif isempty(header.spacing)
        ImageJ_formatted_TIFF.WriteTifStack(output_stack, filepath_output, header.resolution);
    else
        ImageJ_formatted_TIFF.WriteTifStack(output_stack, filepath_output, header.resolution, header.spacing, 'z');
    end
else
    if slices == 1
        if isempty(header.resolution)
            ImageJ_formatted_TIFF.WriteTifStack(output_stack, filepath_output, 0.05, 0.2, 'c');
        else
            ImageJ_formatted_TIFF.WriteTifStack(output_stack, filepath_output, header.resolution, 0.2, 'c');
        end
    else
        output_stack = permute(output_stack, [1 2 4 3]);    % Swap channels and slices
        if isempty(header.resolution) || isempty(header.spacing)
            ImageJ_formatted_TIFF.WriteTifStack(output_stack, filepath_output, 0.05, 0.2, 'cz');
        else
            ImageJ_formatted_TIFF.WriteTifStack(output_stack, filepath_output, header.resolution, header.spacing, 'cz');
        end
    end

end

end