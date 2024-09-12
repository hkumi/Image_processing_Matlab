% ====================================================================== %
% Figures of Merit: image analysis
%   - Contrast-to-Noise Ratio (CNR)
%   - Edge Spread Function (ESF) and Line Spread Function (LSF)
%   - Spatial resolution
%   - Modulation Transfer Function (MTF)
% ====================================================================== %

clc;
clear;
warning('off', 'all');

%% Load Image
[filename, path] = uigetfile('*.tif*', 'Select TIF file', 'MultiSelect', 'off');
if isequal(filename, 0)
    error('No file selected. Exiting...');
end
file = fullfile(path, filename);

image = im2double(imread(file));  % Convert to double for precision
image = imadjust(image);          % Adjust contrast
figure(1);
imshow(image);
title('Original Image');

%% Get Pixel-to-mm Conversion
% Option to provide or calculate pixel-to-mm conversion dynamically
choice = questdlg('Do you want to calculate the px-to-mm ratio from an object in the image?', ...
    'Pixel-to-mm Conversion', 'Yes', 'No, use default', 'No, use default');
if strcmp(choice, 'Yes')
    px2mm = getScaleFromObject(image); % Function to dynamically calculate px2mm
else
    px2mm = 0.08989; % Default conversion factor (you may need to adjust this)
end
disp(['Pixel-to-mm conversion factor: ', num2str(px2mm), ' mm/px']);

%% Image Processing (Optional)
% Dynamically select whether to apply image processing (e.g., filtering)
choice = questdlg('Do you want to apply image processing?', 'Processing Options', 'Yes', 'No', 'No');
if strcmp(choice, 'Yes')
    image = medfilt2(image);      % Median filter to reduce noise
    image = imsharpen(image);     % Sharpen edges
    disp('Image processed with median filter and sharpening.');
end

%% Contrast-to-Noise Ratio (CNR) Calculation
image = imcomplement(image);  % Negative of image for better contrast analysis

[CNR, SNR, Ns] = calculateCNR(image);  % Call function to calculate CNR

image = imcomplement(image);  % Revert image

%% Edge Spread Function (ESF) and Spatial Resolution Calculation
[LSF, FWHM] = calculateESFandLSF(image);

%% Modulation Transfer Function (MTF) Calculation
[MTF, fX, SR_MTF50, SR_MTF10] = calculateMTF(LSF, px2mm);

%% Display Results
disp('------------------------');
disp('Contrast parameters');
disp(['  Contrast-to-Noise Ratio (CNR): ', num2str(CNR)]);
disp(['  Signal-to-Noise Ratio (SNR): ', num2str(SNR)]);
disp(['  Noise dispersion at signal (%): ', num2str(Ns)]);
disp(' ');
disp('Resolution parameters');
disp(['  Spatial resolution (FWHM, px): ', num2str(FWHM)]);
disp(['  Spatial resolution (FWHM, mm): ', num2str(FWHM * px2mm)]);
disp(['  Spatial resolution (MTF@50, mm): ', num2str(SR_MTF50)]);
disp(['  Spatial resolution (MTF@10, mm): ', num2str(SR_MTF10)]);
disp('------------------------');

%% Save Results to File
saveResults(CNR, SNR, Ns, FWHM, SR_MTF50, SR_MTF10);  % Call function to save results

%% Supporting Functions

% Function to calculate Contrast-to-Noise Ratio (CNR)
function [CNR, SNR, Ns] = calculateCNR(image)
    figure(2);
    imshow(image);
    title('Select ROI for Signal');
    
    % Define ROI for signal
    sROI = drawrectangle();
    sM = createMask(sROI);
    avS = mean(image(sM));
    stdS = std(image(sM));

    title('Select ROI for Background');
    % Define ROI for background
    bROI = drawrectangle();
    bM = createMask(bROI);
    avB = mean(image(bM));
    stdB = std(image(bM));

    % Calculate CNR, SNR, and Noise
    CNR = abs(avS - avB) / stdB;
    SNR = avS / stdS;
    Ns = stdS / avS * 100;
end

% Function to calculate Edge Spread Function (ESF) and Line Spread Function (LSF)
function [LSF, FWHM] = calculateESFandLSF(image)
    figure(1);
    imshow(image);
    title('Select ROI for Edge Spread Function');

    % Define ROI for edge
    eROI = drawrectangle('Color', 'red');
    eR = imcrop(image, eROI.Position);

    % Average intensity along the vertical or horizontal direction
    if eROI.AspectRatio > 1
        edge = mean(eR, 2);  % Vertical ROI (horizontal edge)
    else
        edge = mean(eR, 1);  % Horizontal ROI (vertical edge)
    end

    % Ensure intensity profile is low -> high
    if edge(1) > edge(end)
        edge = flip(edge);
    end

    % Fit to error function (Edge Spread Function)
    [xData, yData] = prepareCurveData([], edge);
    ft = fittype('a+b*erf((x-c)/d)', 'independent', 'x', 'dependent', 'y');
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    [ESF_fitr, ~] = fit(xData, yData, ft, opts);

    % Calculate LSF from ESF and normalize
    LSF = differentiate(ESF_fitr, xData);
    LSF = LSF / max(LSF);

    figure(3);
    plot(xData, yData, 'b', 'LineWidth', 1.5);
    hold on;
    plot(ESF_fitr);
    plot(LSF);
    hold off;
    xlabel('Size (px)', 'FontSize', 18);
    ylabel('Intensity');
    legend('Edge Profile', 'ESF Fit', 'LSF');
    title('Edge Spread Function and Line Spread Function');

    % Calculate spatial resolution (Full Width at Half Maximum - FWHM)
    [xData, yData] = prepareCurveData([], LSF);
    ft = fittype('gauss1');
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    LSF_fitr = fit(xData, yData, ft, opts);
    sigma = LSF_fitr.c1;
    FWHM = 2 * sqrt(2 * log(2)) * sigma;  % Full width at half maximum
end

% Function to calculate Modulation Transfer Function (MTF)
function [MTF, fX, SR_MTF50, SR_MTF10] = calculateMTF(LSF, px2mm)
    MTF = abs(fft(LSF));
    MTF = MTF / max(MTF);  % Normalize

    dx = px2mm;            % Spatial step
    df = 1 / dx;           % Spatial frequency step

    % Spatial frequency range
    fX = df * (0:length(LSF)-1) / length(LSF);
    fX(fX > df / 2) = fX(fX > df / 2) - df;

    figure(4);
    plot(fX(1:15), MTF(1:15), 'LineWidth', 2);
    hold on;
    yline(0.5, 'r--', 'MTF 50%');
    yline(0.1, 'b--', 'MTF 10%');
    xlabel('Spatial Frequency (cy/mm)', 'FontSize', 14);
    ylabel('MTF');
    title('Modulation Transfer Function');
    grid on;
    hold off;

    % Calculate spatial frequencies at 50% and 10% MTF
    xMTF_50 = interp1(MTF(1:15), fX(1:15), 0.5);
    xMTF_10 = interp1(MTF(1:15), fX(1:15), 0.1);
    SR_MTF50 = 1 / (2 * xMTF_50);
    SR_MTF10 = 1 / (2 * xMTF_10);
end

% Function to dynamically calculate pixel-to-mm conversion
function px2mm = getScaleFromObject(image)
    figure(1);
    imshow(image);
    title('Select a reference object in the image');

    % Define ROI around reference object
    refROI = drawrectangle();
    pixelLength = refROI.Position(3);  % Width in pixels

    % Ask the user for the real-world size in mm
    prompt = {'Enter the real-world size of the object (mm):'};
    dlgtitle = 'Input';
    dims = [1 35];
    realSize = inputdlg(prompt, dlgtitle, dims);
    realSize = str2double(realSize{1});

    % Calculate px-to-mm conversion
    px2mm = realSize / pixelLength;
end

% Function to save results to a file
function saveResults(CNR, SNR, Ns, FWHM, SR_MTF50, SR_MTF10)
    % Create table of results
    results = table(CNR, SNR, Ns, FWHM, SR_MTF50, SR_MTF10, ...
        'VariableNames', {'CNR', 'SNR', 'Noise', 'FWHM_px', 'MTF50_mm', 'MTF10_mm'});
    
    % Prompt the user to save results
    [file, path] = uiputfile('results.csv', 'Save Results As');
    if isequal(file, 0)
        disp('Save operation canceled.');
    else
        writetable(results, fullfile(path, file));
        disp(['Results saved to: ', fullfile(path, file)]);
    end
end
