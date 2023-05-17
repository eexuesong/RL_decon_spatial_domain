clc;
clear all;
close all;


%% GUI part
[filename_psf, path_psf] = uigetfile('*.tif', 'Choose one PSF tiff file.');
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
ImageJ_formatted_TIFF.WriteTifStack(image_Max_proj, 'image_max_projection.tif', resolution, spacing);
image_bw = single(imbinarize(image_Max_proj / max(image_Max_proj(:)), 'global'));   % Binarization: Calculate global image threshold using Otsu's method
ImageJ_formatted_TIFF.WriteTifStack(image_bw, 'image_bw.tif', resolution, spacing);

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
    singleBead = single(image(xyPos(1) - range(1) : xyPos(1) + range(1), xyPos(2) - range(2) : xyPos(2) + range(2), :));
    [FWHM_x, FWHM_y, FWHM_z] = fwhm_PSF(singleBead, resolution, spacing, 0, 0);
    FWHM_x_nofitting(counter, 1) = FWHM_x;
    FWHM_y_nofitting(counter, 1) = FWHM_y;
    FWHM_z_nofitting(counter, 1) = FWHM_z;
    [FWHM_x, FWHM_y, FWHM_z] = fwhm_PSF(singleBead, resolution, spacing, 0, 1);
    FWHM_x_spinefitting(counter, 1) = FWHM_x;
    FWHM_y_spinefitting(counter, 1) = FWHM_y;
    FWHM_z_spinefitting(counter, 1) = FWHM_z;
    [FWHM_x, FWHM_y, FWHM_z] = fwhm_PSF(singleBead, resolution, spacing, 0, 2);
    FWHM_x_gaussianfitting(counter, 1) = FWHM_x;
    FWHM_y_gaussianfitting(counter, 1) = FWHM_y;
    FWHM_z_gaussianfitting(counter, 1) = FWHM_z;
end
    
disp(['Total successful PSFs: ', num2str(counter)]);


