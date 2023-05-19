%% Function: calculate FWHM of PSF
function [FWHM_x, FWHM_y, FWHM_z] = fwhm_PSF(PSF, pixelSize, pixelDepth, cFlag, fitFlag)
% Feed back the full width at half maximun of the input PSF
% fwhm.m and mygaussfit.m are needed
% cFlag
%       0: use maximum's position as PSF center position
%       1: use matrix's center position as PSF center position
% fitFlag
%       0: no fitting before calculate FWHM
%       1: spine fitting before calculate FWHM
%       2: gaussian fitting before calculate FWHM
%
if(nargin == 1)
    pixelSize = 1;
    pixelDepth = 1;
    cFlag = 0;
    fitFlag = 0;
end

if(nargin == 2)
    pixelDepth = 1;
    cFlag = 0;
    fitFlag = 0;
end

if(nargin == 3)
    cFlag = 0;
    fitFlag = 0;
end

if(nargin == 4)
    fitFlag = 0;
end

% PSF = PSF - mean(PSF(:));
PSF = PSF / max(PSF, [], 'all');
[ny, nx, nz] = size(PSF);
if (ny == 1) || (nx == 1)
    % 1D input
    x = 1:max(ny, nx);
    x = x';
    y = PSF(:);
    if fitFlag == 1
        % fitFlag = 1: spine fitting before calculate FWHM
        xq = 1:0.1:nx;
        yq = interp1(x, y, xq, 'spline');
        FWHM_x = fwhm(xq, yq);
    elseif fitFlag == 2
        % fitFlag = 2: gaussian fitting before calculate FWHM
        [sig, ~] = mygaussfit(x, y);
        %         FWHM_x = sig * 2.3548;
        FWHM_x = sig * 2 * sqrt(2 * log(2));
    else
        % fitFlag = 0: no fitting before calculate FWHM
        FWHM_x = fwhm(x, y);
    end

    FWHM_y = 0;
    FWHM_z = 0;
elseif nz == 1
    % 2D input
    if(cFlag)
        % cFlag = 1: use matrix's center position as PSF center position
        indy = floor((ny + 1) / 2);
        indx = floor((nx + 1) / 2);
    else
        % cFlag = 0: use maximum's position as PSF center position
        [~, ind] = max(PSF(:)); % find maximum value and position
        [indy, indx] = ind2sub([ny, nx], ind(1));
    end

    % Y dimension
    x = 1:ny;
    x = x';
    y = PSF(:, indx);
    y = y(:);
    if fitFlag == 1
        xq = 1:0.1:ny;
        yq = interp1(x, y, xq, 'spline');
        FWHM_y = fwhm(xq, yq);
    elseif fitFlag == 2
        [sig, ~] = mygaussfit(x, y);
        %         FWHM_y = sig * 2.3548;
        FWHM_y = sig * 2 * sqrt(2 * log(2));
    else
        FWHM_y = fwhm(x, y);
    end

    % X dimension
    x = 1:nx;
    x = x';
    y = PSF(indy, :);
    y = y(:);
    if fitFlag == 1
        xq = 1:0.1:nx;
        yq = interp1(x, y, xq, 'spline');
        FWHM_x = fwhm(xq, yq);
    elseif fitFlag == 2
        [sig, ~] = mygaussfit(x, y);
        %         FWHM_x = sig * 2.3548;
        FWHM_x = sig * 2 * sqrt(2 * log(2));
    else
        FWHM_x = fwhm(x, y);
    end

    FWHM_z = 0;
else
    % 3D input
    if(cFlag)
        indy = floor((ny + 1) / 2);
        indx = floor((nx + 1) / 2);
        indz = floor((nz + 1) / 2);
    else
        [~, ind] = max(PSF(:)); % find maximum value and position
        [indy, indx, indz] = ind2sub([ny, nx, nz], ind(1));
    end

    % Y dimension
    x = 1:ny;
    x = x';
    y = PSF(:, indx, indz);
    y = y(:);
    if fitFlag == 1
        xq = 1:0.1:ny;
        yq = interp1(x, y, xq, 'spline');
        FWHM_y = fwhm(xq, yq);
    elseif fitFlag == 2
        [sig, ~] = mygaussfit(x, y);
        %         FWHM_y = sig * 2.3548;
        FWHM_y = sig * 2 * sqrt(2 * log(2));
    else
        FWHM_y = fwhm(x, y);
    end

    % X dimension
    x = 1:nx;
    x = x';
    y = PSF(indy, :, indz);
    y = y(:);
    if fitFlag == 1
        xq = 1:0.1:nx;
        yq = interp1(x, y, xq, 'spline');
        FWHM_x = fwhm(xq, yq);
    elseif fitFlag == 2
        [sig, ~] = mygaussfit(x, y);
        %         FWHM_x = sig * 2.3548;
        FWHM_x = sig * 2 * sqrt(2 * log(2));
    else
        FWHM_x = fwhm(x, y);
    end

    % Z dimension
    x = 1:nz;
    x = x';
    y = PSF(indy, indx, :);
    y = y(:);
    if fitFlag == 1
        xq = 1:0.1:nz;
        yq = interp1(x, y, xq, 'spline');
        FWHM_z = fwhm(xq, yq);
    elseif fitFlag == 2
        [sig, ~] = mygaussfit(x, y);
        %         FWHM_z = sig * 2.3548;
        FWHM_z = sig * 2 * sqrt(2 * log(2));
    else
        FWHM_z = fwhm(x, y);
    end
end

FWHM_x = FWHM_x * pixelSize;
FWHM_y = FWHM_y * pixelSize;
FWHM_z = FWHM_z * pixelDepth;
end