% Xuesong Li 05/10/2023:

%% GUI part
% load raw data
[filename_data, path_data] = uigetfile('*.tif', 'Choose any one of raw data');

% set parameters
inputDialog = uifigure('Name', 'Set Parameters', 'Position', [300 500 500 320]);

% panels
pnl_60x = uipanel(inputDialog, 'Position', [10 210 480 105]);
pnl_100x = uipanel(inputDialog, 'Position', [10 100 480 105]);

% set FWHM of PSFs (60X)
uilabel(pnl_60x, 'Text', 'UPLSAPO60XS2', 'Position', [10 80 100 20]);
uilabel(pnl_60x, 'Text', '405 nm:', 'Position', [90 60 60 20]);
uilabel(pnl_60x, 'Text', '488 nm:', 'Position', [190 60 60 20]);
uilabel(pnl_60x, 'Text', '561 nm:', 'Position', [290 60 60 20]);
uilabel(pnl_60x, 'Text', '633 nm:', 'Position', [390 60 60 20]);

uilabel(pnl_60x, 'Text', 'XY (nm):', 'Position', [30 40 60 20]);
txa_60x_FWHM_xy{1} = uitextarea(pnl_60x, 'Value', 'NA', 'Position', [90 40 60 20]);
txa_60x_FWHM_xy{2} = uitextarea(pnl_60x, 'Value', 'NA', 'Position', [190 40 60 20]);
txa_60x_FWHM_xy{3} = uitextarea(pnl_60x, 'Value', 'NA', 'Position', [290 40 60 20]);
txa_60x_FWHM_xy{4} = uitextarea(pnl_60x, 'Value', 'NA', 'Position', [390 40 60 20]);

uilabel(pnl_60x, 'Text', 'Z (nm):', 'Position', [38 10 60 20]);
txa_60x_FWHM_z{1} = uitextarea(pnl_60x, 'Value', 'NA', 'Position', [90 10 60 20]);
txa_60x_FWHM_z{2} = uitextarea(pnl_60x, 'Value', 'NA', 'Position', [190 10 60 20]);
txa_60x_FWHM_z{3} = uitextarea(pnl_60x, 'Value', 'NA', 'Position', [290 10 60 20]);
txa_60x_FWHM_z{4} = uitextarea(pnl_60x, 'Value', 'NA', 'Position', [390 10 60 20]);

% set FWHM of PSFs (100X)
uilabel(pnl_100x, 'Text', 'UPLSAPO100XO', 'Position', [10 80 100 20]);
uilabel(pnl_100x, 'Text', '405 nm:', 'Position', [90 60 60 20]);
uilabel(pnl_100x, 'Text', '488 nm:', 'Position', [190 60 60 20]);
uilabel(pnl_100x, 'Text', '561 nm:', 'Position', [290 60 60 20]);
uilabel(pnl_100x, 'Text', '633 nm:', 'Position', [390 60 60 20]);

uilabel(pnl_100x, 'Text', 'XY (nm):', 'Position', [30 40 60 20]);
txa_100x_FWHM_xy{1} = uitextarea(pnl_100x, 'Value', '155', 'Position', [90 40 60 20]);
txa_100x_FWHM_xy{2} = uitextarea(pnl_100x, 'Value', '184', 'Position', [190 40 60 20]);
txa_100x_FWHM_xy{3} = uitextarea(pnl_100x, 'Value', '207', 'Position', [290 40 60 20]);
txa_100x_FWHM_xy{4} = uitextarea(pnl_100x, 'Value', '240', 'Position', [390 40 60 20]);

uilabel(pnl_100x, 'Text', 'Z (nm):', 'Position', [38 10 60 20]);
txa_100x_FWHM_z{1} = uitextarea(pnl_100x, 'Value', '350', 'Position', [90 10 60 20]);
txa_100x_FWHM_z{2} = uitextarea(pnl_100x, 'Value', '400', 'Position', [190 10 60 20]);
txa_100x_FWHM_z{3} = uitextarea(pnl_100x, 'Value', '450', 'Position', [290 10 60 20]);
txa_100x_FWHM_z{4} = uitextarea(pnl_100x, 'Value', '533', 'Position', [390 10 60 20]);

% set iteration number
uilabel(inputDialog, 'Text', 'Iter Num:', 'Position', [40 60 60 20]);
txa_itNum{1} = uitextarea(inputDialog, 'Value', '10', 'Position', [100 60 60 20]);
txa_itNum{2} = uitextarea(inputDialog, 'Value', '10', 'Position', [200 60 60 20]);
txa_itNum{3} = uitextarea(inputDialog, 'Value', '10', 'Position', [300 60 60 20]);
txa_itNum{4} = uitextarea(inputDialog, 'Value', '10', 'Position', [400 60 60 20]);

% buttons to confirm/cancel processing
pb1 = uibutton(inputDialog, 'push', 'Position', [180 20 50 20], 'Text', 'Yes',...
    'ButtonPushedFcn', @(pb1, event) pbYes(inputDialog, txa_60x_FWHM_xy, txa_60x_FWHM_z, txa_100x_FWHM_xy, txa_100x_FWHM_z, txa_itNum, path_data));
pb2 = uibutton(inputDialog, 'push', 'Position', [280 20 50 20], 'Text', 'Cancel',...
    'ButtonPushedFcn', @(pb1, event) pbCancel(inputDialog));

%% Cancel dialog box
function pbCancel(inputDialog)
    delete(inputDialog);
    disp('Processing cancelled!!!')
end

%% Yes dialog box
function pbYes(inputDialog, txa_60x_FWHM_xy, txa_60x_FWHM_z, txa_100x_FWHM_xy, txa_100x_FWHM_z, txa_itNum, path_data)
%     set(inputDialog, 'Pointer', 'watch');
    myFiles = dir(path_data);
    filenames = {myFiles.name};
    mask = endsWith(filenames, {'.tif', '.tiff'}, 'IgnoreCase', true);
    ListOfImages = filenames(mask);
    
    for file_index = 1:length(ListOfImages)
        filename = strcat(path_data, ListOfImages{file_index});
        header = ImageJ_formatted_TIFF.get_header(filename);

        if ~isempty(header.IJMetadata)
            FWHM_x = zeros([1, length(header.IJMetadata.ChNames)]);
            FWHM_y = zeros([1, length(header.IJMetadata.ChNames)]);
            FWHM_z = zeros([1, length(header.IJMetadata.ChNames)]);
            itNum = zeros([1, length(header.IJMetadata.ChNames)]);

            for channel_index = 1:length(header.IJMetadata.ChNames)
                if contains(header.MicroManagerMetadata, 'UPLSAPO60XS2')
                    switch header.IJMetadata.ChNames{channel_index}
                        case 'GFP'
                            FWHM_x(channel_index) = str2double(get(txa_60x_FWHM_xy{2}, 'Value'));
                            FWHM_y(channel_index) = str2double(get(txa_60x_FWHM_xy{2}, 'Value'));
                            FWHM_z(channel_index) = str2double(get(txa_60x_FWHM_z{2}, 'Value'));
                            itNum(channel_index) = str2double(get(txa_itNum{2}, 'Value'));
                        case 'mCherry'
                            FWHM_x(channel_index) = str2double(get(txa_60x_FWHM_xy{3}, 'Value'));
                            FWHM_y(channel_index) = str2double(get(txa_60x_FWHM_xy{3}, 'Value'));
                            FWHM_z(channel_index) = str2double(get(txa_60x_FWHM_z{3}, 'Value'));
                            itNum(channel_index) = str2double(get(txa_itNum{3}, 'Value'));
                    end
                elseif contains(header.MicroManagerMetadata, 'UPLSAPO100XO')
                    switch header.IJMetadata.ChNames{channel_index}
                        case 'GFP'
                            FWHM_x(channel_index) = str2double(get(txa_100x_FWHM_xy{2}, 'Value'));
                            FWHM_y(channel_index) = str2double(get(txa_100x_FWHM_xy{2}, 'Value'));
                            FWHM_z(channel_index) = str2double(get(txa_100x_FWHM_z{2}, 'Value'));
                            itNum(channel_index) = str2double(get(txa_itNum{2}, 'Value'));
                        case 'mCherry'
                            FWHM_x(channel_index) = str2double(get(txa_100x_FWHM_xy{3}, 'Value'));
                            FWHM_y(channel_index) = str2double(get(txa_100x_FWHM_xy{3}, 'Value'));
                            FWHM_z(channel_index) = str2double(get(txa_100x_FWHM_z{3}, 'Value'));
                            itNum(channel_index) = str2double(get(txa_itNum{3}, 'Value'));
                    end
                else
                    error('iSIM_decon_ui: no objective information is found');
                end
                
                xres = FWHM_x / 1000 / double(header.resolution);
                yres = FWHM_y / 1000 / double(header.resolution);
                zres = FWHM_z / 1000 / abs(header.spacing);
            end

            gpu_decon_3d(xres, yres, zres, itNum, filename);
        end  
    end

    delete(inputDialog);
end