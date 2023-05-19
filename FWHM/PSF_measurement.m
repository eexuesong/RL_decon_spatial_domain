clc;
clear all;
close all;

verbose = 1;


%% GUI part
[filename_psf, path_psf] = uigetfile('*.tif', 'Choose one PSF tiff file.');

%% Create an output folder
path_output = strcat(path_psf, '\', 'results\');
mkdir(path_output);
if verbose
    path_no_fitting = strcat(path_output, 'no_fitting\');
    path_spine = strcat(path_output, 'spine_interpolation\');
    path_gaussian = strcat(path_output, 'Gaussian_fitting\');
    mkdir(path_no_fitting);
    mkdir(path_spine);
    mkdir(path_gaussian);
end

%% Read an PSF image
[image, header] = ImageJ_formatted_TIFF.ReadTifStack(fullfile(path_psf, filename_psf));
image = single(image);

if isempty(header.resolution)
    resolution = 0.05;             % unit: um
else
    resolution = header.resolution;
end

if isempty(header.spacing)
    spacing = 0.125;                 % unit: um
else
    spacing = header.spacing;
end

background = 100;

%% Background subtraction
image = image - background;
image(image < 0) = 0;

%% Find beads
image_Max_proj = max(image, [], 3);     % Maxinum projection
ImageJ_formatted_TIFF.WriteTifStack(image_Max_proj, fullfile(path_output, 'image_max_projection.tif'), resolution, spacing);
image_bw = single(imbinarize(image_Max_proj / max(image_Max_proj(:)), 'global'));   % Binarization: Calculate global image threshold using Otsu's method
ImageJ_formatted_TIFF.WriteTifStack(image_bw, fullfile(path_output, 'image_bw.tif'), resolution, spacing);

[l, m] = bwlabel(image_bw);             % Label connected components in 2-D binary image
status = regionprops(l, 'BoundingBox'); % Measure properties of image regions: smallest rectangle containing the region
centroid = regionprops(l, 'Centroid');  % Center of mass of the region
area = regionprops(l, 'Area');

centroids = cat(1, centroid.Centroid);
imshow(image_bw);
hold on;
plot(centroids(:, 1), centroids(:, 2), 'b*')    % centroids(:, 1): x; centroids(:, 2): y
hold off;
mask = double(image_bw);

range = [25, 25, 20];                   % cut a region [y, x, z] [height / 2, width / 2, depth / 2] containing beads
counter = 0;

All_coordinates = zeros(length(centroid), 2);
for i = 1 : length(centroid)
   mGlobalCoor = round(centroid(i).Centroid(1, 2:-1:1));    % Swap the (x, y) coordinates to match Matlab order (y, x)
   All_coordinates(i, :) = mGlobalCoor;
end

Diff_coordinates = zeros(length(centroid), 2);
for i = 1 : length(centroid)
    % Bypass those areas which are too small or too large
    if (area(i).Area(1) < 5) || (area(i).Area(1) > 200)  
        continue;
    end
    
    mGlobalCoor = round(centroid(i).Centroid(1, 2:-1:1));
    
    % Make sure that all the other regions are far enough away from the current region
    Diff_coordinates = All_coordinates - All_coordinates(i, :);
    Diff_coordinates(i, :) = [];        % Remove the current i-th row;
    Distance = sqrt(Diff_coordinates(:, 1).^2 + Diff_coordinates(:, 2).^2);
    if min(Distance) < sqrt(2) * range(1)
        continue;
    end
    
    posErr = [2, 2];
    % Make sure that the current region are not close to the edges of the image
    if (mGlobalCoor(1) < posErr(1) + 1) || (mGlobalCoor(2) < posErr(2) + 1) ...
            || (mGlobalCoor(1) > header.ImageLength - posErr(1)) || (mGlobalCoor(2) > header.ImageWidth - posErr(2))
        continue;
    end
    
    % find the maximum intensity in the surrounding region
    localRegion = image_Max_proj(mGlobalCoor(1) - posErr(1) : mGlobalCoor(1) + posErr(1), ...
        mGlobalCoor(2) - posErr(2) : mGlobalCoor(2) + posErr(2));
    [tempy, tempx] =  find(localRegion == max(max(localRegion)));
    xyPos = mGlobalCoor - posErr + [tempy, tempx] - 1;
    z = image(xyPos(1, 1), xyPos(1, 2), :);
    [~, zPos] = max(z);
    
    % Make sure that we can crop the current PSF with enough pixels
    if (xyPos(1) < range(1) + 1) || (xyPos(1) > header.ImageLength - range(1)) || ...
            (xyPos(2) < range(2) + 1) || (xyPos(2) > header.ImageLength - range(2)) || ...
            (zPos < range(3) + 1) || (zPos > header.images - range(3))
        continue;
    else
        counter = counter + 1;
    end
    
    % Crop the image stack of single bead
    singleBead = single(image(xyPos(1) - range(1) : xyPos(1) + range(1), xyPos(2) - range(2) : xyPos(2) + range(2), zPos - range(3) : zPos + range(3)));
    singleBead = singleBead / max(singleBead, [], 'all');
    [ny_singleBead, nx_singleBead, nz_singleBead] = size(singleBead);

    % FWHM measurement
    [FWHM_x, FWHM_y, FWHM_z] = fwhm_PSF(singleBead, resolution, spacing, 0, 0);
    FWHM_x_nofitting(counter, 1) = FWHM_x;
    FWHM_y_nofitting(counter, 1) = FWHM_y;
    FWHM_z_nofitting(counter, 1) = FWHM_z;
    [FWHM_x, FWHM_y, FWHM_z] = fwhm_PSF(singleBead, resolution, spacing, 0, 1);
    FWHM_x_spine(counter, 1) = FWHM_x;
    FWHM_y_spine(counter, 1) = FWHM_y;
    FWHM_z_spine(counter, 1) = FWHM_z;
    [FWHM_x, FWHM_y, FWHM_z] = fwhm_PSF(singleBead, resolution, spacing, 0, 2);
    FWHM_x_gaussian(counter, 1) = FWHM_x;
    FWHM_y_gaussian(counter, 1) = FWHM_y;
    FWHM_z_gaussian(counter, 1) = FWHM_z;
    
    if verbose
        %% Save each PSF image
        [~, ind] = max(singleBead(:)); % find maximum value and position
        [indy, indx, indz] = ind2sub(size(singleBead), ind(1));
        singleBead_xy = singleBead(:, :, indz);
        ImageJ_formatted_TIFF.WriteTifStack(singleBead_xy, fullfile(path_output, strcat('psf_', num2str(counter), '_xy.tif')), resolution, spacing);
        singleBead_xz = squeeze(singleBead(indy, :, :));
        singleBead_xz = imresize(singleBead_xz, [nx_singleBead, round(nz_singleBead * spacing / resolution)]);
        singleBead_xz = permute(singleBead_xz, [2 ,1]);
        ImageJ_formatted_TIFF.WriteTifStack(singleBead_xz, fullfile(path_output, strcat('psf_', num2str(counter), '_xz.tif')), resolution, spacing);

        %% Save Y profile of each PSF
        x = (-range(1):range(1))' * resolution;
        y = singleBead(:, indx, indz);
        y = y(:);

        % no fitting
        [rail, lead] = fwhm_2(x, y);
        fig = figure('visible','off');
        hold on;
        plot(x, y, '-ok');
        plot([lead, lead], [0, 0.5], '--r');
        plot([rail, rail], [0, 0.5], '--r');
        text(rail, 0.5, strcat('\leftarrow FWHM (no fitting) = ', num2str(rail - lead), ' \mum'), 'Color', 'r');
        ylim([0, 1.05]);
        xlabel('Y distance (\mum)');
        ylabel('Intensity (a.u.)');
        hold off;
        saveas(fig, fullfile(path_no_fitting, strcat('FWHM_', num2str(counter), '_y.tif')));
        close(fig);

        % Spine interpolation
        xq = (-range(1):0.1:range(1))' * resolution;
        yq = interp1(x, y, xq, 'spline');
        [rail, lead] = fwhm_2(xq, yq);
        fig = figure('visible','off');
        hold on;
        plot(x, y, 'ok');
        plot(xq, yq, '-k');
        plot([lead, lead], [0, 0.5], '--r');
        plot([rail, rail], [0, 0.5], '--r');
        text(rail, 0.5, strcat('\leftarrow FWHM (spine interpolation) = ', num2str(rail - lead), ' \mum'), 'Color', 'r');
        ylim([0, 1.05]);
        xlabel('Y distance (\mum)');
        ylabel('Intensity (a.u.)');
        hold off;
        saveas(fig, fullfile(path_spine, strcat('FWHM_', num2str(counter), '_y.tif')));
        close(fig);

        % Gaussian fitting
        [sigma, mu, A] = mygaussfit(x, y);
        x_gaussian = (-range(1):0.1:range(1))' * resolution;
        y_gaussian = A * exp(-(x_gaussian - mu).^2 / (2 * sigma^2));
        fwhm_gaussian = sigma * 2 * sqrt(2 * log(2));
        fig = figure('visible','off');
        hold on;
        plot(x, y, 'ok');
        plot(x_gaussian, y_gaussian, '-r');
        plot([mu - fwhm_gaussian / 2, mu - fwhm_gaussian / 2], [0, A / 2], '--r');
        plot([mu + fwhm_gaussian / 2, mu + fwhm_gaussian / 2], [0, A / 2], '--r');
        text(mu + fwhm_gaussian / 2, A / 2, strcat('\leftarrow FWHM (Gaussian fitting) = ', num2str(fwhm_gaussian), ' \mum'), 'Color', 'r');
        ylim([0, 1.05]);
        xlabel('Y distance (\mum)');
        ylabel('Intensity (a.u.)');
        hold off;
        saveas(fig, fullfile(path_gaussian, strcat('FWHM_', num2str(counter), '_y.tif')));
        close(fig);

        %% Save X profile of each PSF
        x = (-range(2):range(2))' * resolution;
        y = singleBead(indy, :, indz);
        y = y(:);

        % no fitting
        [rail, lead] = fwhm_2(x, y);
        fig = figure('visible','off');
        hold on;
        plot(x, y, '-ok');
        plot([lead, lead], [0, 0.5], '--r');
        plot([rail, rail], [0, 0.5], '--r');
        text(rail, 0.5, strcat('\leftarrow FWHM (no fitting) = ', num2str(rail - lead), ' \mum'), 'Color', 'r');
        ylim([0, 1.05]);
        xlabel('X distance (\mum)');
        ylabel('Intensity (a.u.)');
        hold off;
        saveas(fig, fullfile(path_no_fitting, strcat('FWHM_', num2str(counter), '_x.tif')));
        close(fig);

        % Spine interpolation
        xq = (-range(2):0.1:range(2))' * resolution;
        yq = interp1(x, y, xq, 'spline');
        [rail, lead] = fwhm_2(xq, yq);
        fig = figure('visible','off');
        hold on;
        plot(x, y, 'ok');
        plot(xq, yq, '-k');
        plot([lead, lead], [0, 0.5], '--r');
        plot([rail, rail], [0, 0.5], '--r');
        text(rail, 0.5, strcat('\leftarrow FWHM (spine interpolation) = ', num2str(rail - lead), ' \mum'), 'Color', 'r');
        ylim([0, 1.05]);
        xlabel('X distance (\mum)');
        ylabel('Intensity (a.u.)');
        hold off;
        saveas(fig, fullfile(path_spine, strcat('FWHM_', num2str(counter), '_x.tif')));
        close(fig);

        % Gaussian fitting
        [sigma, mu, A] = mygaussfit(x, y);
        x_gaussian = (-range(2):0.1:range(2))' * resolution;
        y_gaussian = A * exp(-(x_gaussian - mu).^2 / (2 * sigma^2));
        fwhm_gaussian = sigma * 2 * sqrt(2 * log(2));
        fig = figure('visible','off');
        hold on;
        plot(x, y, 'ok');
        plot(x_gaussian, y_gaussian, '-r');
        plot([mu - fwhm_gaussian / 2, mu - fwhm_gaussian / 2], [0, A / 2], '--r');
        plot([mu + fwhm_gaussian / 2, mu + fwhm_gaussian / 2], [0, A / 2], '--r');
        text(mu + fwhm_gaussian / 2, A / 2, strcat('\leftarrow FWHM (Gaussian fitting) = ', num2str(fwhm_gaussian), ' \mum'), 'Color', 'r');
        ylim([0, 1.05]);
        xlabel('X distance (\mum)');
        ylabel('Intensity (a.u.)');
        hold off;
        saveas(fig, fullfile(path_gaussian, strcat('FWHM_', num2str(counter), '_x.tif')));
        close(fig);

        %% Save Z profile of each PSF
        x = (-range(3):range(3))' * spacing;
        y = singleBead(indy, indx, :);
        y = y(:);

        % no fitting
        [rail, lead] = fwhm_2(x, y);
        fig = figure('visible','off');
        hold on;
        plot(x, y, '-ok');
        plot([lead, lead], [0, 0.5], '--r');
        plot([rail, rail], [0, 0.5], '--r');
        text(rail, 0.5, strcat('\leftarrow FWHM (no fitting) = ', num2str(rail - lead), ' \mum'), 'Color', 'r');
        xlabel('Z distance (\mum)');
        ylabel('Intensity (a.u.)');
        hold off;
        saveas(fig, fullfile(path_no_fitting, strcat('FWHM_', num2str(counter), '_z.tif')));
        close(fig);

        % Spine interpolation
        xq = (-range(3):0.1:range(3))' * spacing;
        yq = interp1(x, y, xq, 'spline');
        [rail, lead] = fwhm_2(xq, yq);
        fig = figure('visible','off');
        hold on;
        plot(x, y, 'ok');
        plot(xq, yq, '-k');
        plot([lead, lead], [0, 0.5], '--r');
        plot([rail, rail], [0, 0.5], '--r');
        text(rail, 0.5, strcat('\leftarrow FWHM (spine interpolation) = ', num2str(rail - lead), ' \mum'), 'Color', 'r');
        ylim([0, 1.05]);
        xlabel('Z distance (\mum)');
        ylabel('Intensity (a.u.)');
        hold off;
        saveas(fig, fullfile(path_spine, strcat('FWHM_', num2str(counter), '_z.tif')));
        close(fig);

        % Gaussian fitting
        [sigma, mu, A] = mygaussfit(x, y);
        x_gaussian = (-range(3):0.1:range(3))' * spacing;
        y_gaussian = A * exp(-(x_gaussian - mu).^2 / (2 * sigma^2));
        fwhm_gaussian = sigma * 2 * sqrt(2 * log(2));
        fig = figure('visible','off');
        hold on;
        plot(x, y, 'ok');
        plot(x_gaussian, y_gaussian, '-r');
        plot([mu - fwhm_gaussian / 2, mu - fwhm_gaussian / 2], [0, A / 2], '--r');
        plot([mu + fwhm_gaussian / 2, mu + fwhm_gaussian / 2], [0, A / 2], '--r');
        text(mu + fwhm_gaussian / 2, A / 2, strcat('\leftarrow FWHM (Gaussian fitting) = ', num2str(fwhm_gaussian), ' \mum'), 'Color', 'r');
        ylim([0, 1.05]);
        xlabel('Z distance (\mum)');
        ylabel('Intensity (a.u.)');
        hold off;
        saveas(fig, fullfile(path_gaussian, strcat('FWHM_', num2str(counter), '_z.tif')));
        close(fig);
    end
    
    writematrix(FWHM_x_nofitting, fullfile(path_output, 'fwhm.xlsx'), 'Sheet', 'X', 'Range', strcat('A1:A', num2str(counter)));
    writematrix(FWHM_x_spine, fullfile(path_output, 'fwhm.xlsx'), 'Sheet', 'X', 'Range', strcat('B1:B', num2str(counter)));
    writematrix(FWHM_x_gaussian, fullfile(path_output, 'fwhm.xlsx'), 'Sheet', 'X', 'Range', strcat('C1:C', num2str(counter)));

    writematrix(FWHM_y_nofitting, fullfile(path_output, 'fwhm.xlsx'), 'Sheet', 'Y', 'Range', strcat('A1:A', num2str(counter)));
    writematrix(FWHM_y_spine, fullfile(path_output, 'fwhm.xlsx'), 'Sheet', 'Y', 'Range', strcat('B1:B', num2str(counter)));
    writematrix(FWHM_y_gaussian, fullfile(path_output, 'fwhm.xlsx'), 'Sheet', 'Y', 'Range', strcat('C1:C', num2str(counter)));

    writematrix(FWHM_z_nofitting, fullfile(path_output, 'fwhm.xlsx'), 'Sheet', 'Z', 'Range', strcat('A1:A', num2str(counter)));
    writematrix(FWHM_z_spine, fullfile(path_output, 'fwhm.xlsx'), 'Sheet', 'Z', 'Range', strcat('B1:B', num2str(counter)));
    writematrix(FWHM_z_gaussian, fullfile(path_output, 'fwhm.xlsx'), 'Sheet', 'Z', 'Range', strcat('C1:C', num2str(counter)));
end
    
disp(['Total successful PSFs: ', num2str(counter)]);


