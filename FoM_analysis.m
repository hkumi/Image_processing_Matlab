% ====================================================================== %
% Figures of Merit: image analysis
%   - Contrast-to-Noise Ratio (CNR)
%   - Edge Spread Function (ESF) and Line Spread Function (LSF)
%   - Spatial resolution
%   - Modulation Transfer Function (MTF)
% ====================================================================== %

clc;
clear; %to clear off the workspace 
warning('off', 'all'); % to disable all the warning. 

%% Load image

% Select path of image
[filename, path] = uigetfile('*.tif*', 'Select TIF files', 'MultiSelect', 'off');
if isequal(filename, 0)
    error('No file selected. Exiting...');
end
file = fullfile(path,filename);

% Read file, convert to double rep. and adjust contrast
image = im2double(imread(file));
image = imadjust(image);
%image = imcomplement(image);

% Plot image
figure(1)
imshow(image);

title('Original Image');
%% Parameters
px2mm = 0.08989; %1/20; %mm/px %this basicaly needs to be adjusted with the right parameter because it affects the overal calculation


improc = 0; % 1=process image (ah, medfilt, sharp)

%% Proccess image (if set)

%% Image Processing (Optional)
% Dynamically select whether to apply image processing (e.g., filtering)
choice = questdlg('Do you want to apply image processing?', 'Processing Options', 'Yes', 'No', 'No');
if strcmp(choice, 'Yes')
    image = medfilt2(image);      % Median filter to reduce noise
    image = imsharpen(image);     % Sharpen edges
    disp('Image processed with median filter and sharpening.');
end

%% Contrast-to-Noise Ratio (CNR)

image = imcomplement(image);    % Negative of image
title('Select ROI for Signal');
%Define ROI for signal
sROI = drawrectangle();
sM = createMask(sROI);


avS = mean(image(sM));
stdS = std(image(sM));

%Define ROI for background
title('Select ROI for Background');
bROI = drawrectangle();
bM = createMask(bROI);

avB = mean(image(bM));
stdB = std(image(bM));

%Calculate CNR, SNR and Noise
CNR = abs(avS - avB)/stdB;
SNR = avS/stdS;
Ns = stdS/(avS)*100;

image = imcomplement(image);

%% Edge Spread Function (ESF) and spatial resolution

%Reset image
figure(1)
imshow(image);

%Define ROI for edge
title('Select ROI for Edge Spread Function');
eROI = drawrectangle("Color","red");
eR = imcrop(image,eROI.Position);

if eROI.AspectRatio > 1 %Vertical ROI (horizontal edge)
    edge = mean(eR,2);
    
elseif eROI.AspectRatio < 1 % Horizontal ROI (vertical edge)
    edge = mean(eR,1);
    
end

if edge(1) > edge(end)  % Flip intensity profile to low->high
    edge = flip(edge);
    
end

%Fit to error function (Edge Spread Function)
[xData, yData] = prepareCurveData( [], edge );

ft = fittype( 'a+b*erf((x-c)/d)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';

[ESF_fitr, ESF_gof] = fit( xData, yData, ft, opts );

%Calculate LSF from ESF and normalize
LSF = differentiate(ESF_fitr,xData);
LSF = LSF/max(LSF);

figure(2)
hold on
plot(xData, yData,'b','LineWidth',1.5);
plot(ESF_fitr)
plot(LSF)
hold off
xlabel('Size (px)','FontSize',18)
ylabel('Intensity')
legend('Edge Profile','ESF','LSF')
title('Edge Spread Function and Line Spread Function');

%Calculate resolution from LSF
[xData, yData] = prepareCurveData( [], LSF );

ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

[LSF_fitr, LSF_gof] = fit( xData, yData, ft, opts );

sigma = LSF_fitr.c1;
FWHM = 2*sqrt(2*log(2))*sigma;

%% Modulation Transfer Function (MTF)

MTF = abs(fft(LSF));
MTF = MTF/max(MTF); %Normalize

%fX = linspace(0,length(xData)/2,numel(MTF)/2)/px2mm;

dx = (xData(2) - xData(1))*px2mm;   % Spatial step
df = 1/dx;                          % Spatial freq. step

fX = df*(0:length(xData)-1)/length(xData);     % Spatial freq. range
fX(fX>df/2) = fX(fX>df/2)-df;

%figure(3)
%plot(fX(1:15),MTF(1:15),'LineWidth', 2)
%xlabel('cy/mm');
%ylabel('MTF');
figure(3);
plot(fX(1:15), MTF(1:15), 'LineWidth', 2);
hold on;
yline(0.5, 'r--', 'MTF 50%');
yline(0.1, 'b--', 'MTF 10%');
xlabel('Spatial Frequency (cy/mm)', 'FontSize', 14);
ylabel('MTF');
title('Modulation Transfer Function');
grid on;
hold off;

% Calculate spatial freq. and resolution at 0.5 and 0.1 of MTF
xMTF_50 = interp1(MTF(1:15),fX(1:15),0.5);
xMTF_10 = interp1(MTF(1:15),fX(1:15),0.1);
SR_MTF50 = 1/(2*xMTF_50);
SR_MTF10 = 1/(2*xMTF_10);

%% Display results

disp('------------------------')
disp('Contrast parameters');
disp(strcat('  Contrast-to-Noise Ratio (CNR): ', num2str(CNR)));
disp(strcat('  Signal-to-Noise Ratio (SNR): ', num2str(SNR)));
disp(strcat('  Noise dispersion at signal (%): ', num2str(Ns)));
disp(' ');
disp('Resolution parameters');
disp(strcat('  Spatial resolution (FWHM, px): ', num2str(FWHM)));
disp(strcat('  Spatial resolution (FWHM, mm): ', num2str(FWHM*px2mm)));
disp(strcat('  Spatial resolution (MTF@50, mm): ', num2str(SR_MTF50)));
disp(strcat('  Spatial resolution (MTF@10, mm): ', num2str(SR_MTF10)));
disp('------------------------')
%% Save Results to File
saveResults(CNR, SNR, Ns, FWHM, SR_MTF50, SR_MTF10);  % Call function to save results


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