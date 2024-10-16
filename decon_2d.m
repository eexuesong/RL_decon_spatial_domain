function decon_2d(type, xres, yres, background, iteration, filename)
% Reads a image file, processes it, saves an output tif as 'filename_decon.tif'.

% 09/06/17: gpu based serial decon, using gpu_decon. (USE THIS ONE)

% The filename should not contain space
% xres yres  In pixels, the full width at half maximum of the Gaussian PSF.

% [4, 4]:
% The actual FWHM in X-Y is around 180 nm. Thefore the value = 0.18 / pixel_size = 0.18 / 0.046 = 3.9 pixels.
% X,Y is normally keep constant which is always around 4 pixels.

% For 60x water 1.2 NA lens FWHM in X-Y is [2.5, 2.5].
% For ORCA_Quest, pixel size 46 nm. For 100x oil lens The X-Y should be 4 pixels without binning. And 2 with 2x2 binning.


%% Read tiff stack file
[image, header] = ImageJ_formatted_TIFF.ReadTifStack(filename);
if isempty(header.channels)
    channels = 1;
else
    channels = header.channels;
    if (channels ~= size(xres, 2)) || (channels ~= size(yres, 2)) || (channels ~= size(background, 2)) || (channels ~= size(iteration, 2))
        error("Channel number of image stack is inconsistent with xres, yres, background or iteration array size.");
    end
end

% 2D slice or 3D stack
if isempty(header.frames)
    frames = 1;
else
    frames = header.frames;
end

output_stack = zeros(size(image), 'uint16');
if (channels ~= 1) && (frames ~= 1)
    % 4D 'ct' stack only
    output_stack = permute(output_stack, [1 2 4 3]); % Swap channels and frames
end


%% RL deconvolution part
for channel_index = 1:channels
    if channels == 1
        input_image = image;
        input_image = input_image - background;
    else
        input_image = image(:, :, channel_index, :);
        input_image = input_image - background(channel_index);
        if frames > 1
            input_image = squeeze(permute(input_image, [1 2 4 3]));    % Swap frames and channels
        end
    end

    input_image(input_image < 0) = 0;

%     if min(input_image, [], 'all') <= 0
%         input_image = input_image + 0.00001;
%     end
       
    zres = 0;   % zres is meaningless here.
    FWHM = [xres(channel_index), yres(channel_index), zres(channel_index)];

    for frame_index = 1:frames
        disp(['Frame index: ', num2str(frame_index)]);
        if type == "gpu"
            output_image = decon_gpu(input_image(:, :, frame_index), FWHM, iteration(channel_index));
        elseif type == "cpu"
            output_image = decon_cpu(input_image(:, :, frame_index), FWHM, iteration(channel_index));
        else
            output_image = decon_core(input_image(:, :, frame_index), FWHM, iteration(channel_index));
        end

        if (header.BitsPerSample == 16)
            if max(output_image, [], 'all') <= 65535
                output_image = uint16(output_image);
            else
                output_image = uint16(65535 * output_image ./ max(output_image, [], 'all'));
            end
        end

        if frames == 1
            output_stack(:, :, channel_index) = output_image;
        else
            output_stack(:, :, frame_index, channel_index) = output_image;
        end
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
        ImageJ_formatted_TIFF.WriteTifStack(output_stack, filepath_output, header.resolution, header.spacing, 't');
    end
else
    if frames == 1
        if isempty(header.resolution)
            ImageJ_formatted_TIFF.WriteTifStack(output_stack, filepath_output, 0.05, 0.2, 'c');
        else
            ImageJ_formatted_TIFF.WriteTifStack(output_stack, filepath_output, header.resolution, 0.2, 'c');
        end
    else
        output_stack = permute(output_stack, [1 2 4 3]);    % Swap channels and frames
        if isempty(header.resolution) || isempty(header.spacing)
            ImageJ_formatted_TIFF.WriteTifStack(output_stack, filepath_output, 0.05, 0.2, 'ct');
        else
            ImageJ_formatted_TIFF.WriteTifStack(output_stack, filepath_output, header.resolution, header.spacing, 'ct');
        end
    end

end

end