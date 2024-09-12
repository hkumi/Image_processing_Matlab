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

%% Load image

% Select path of image
[filename, path] = uigetfile('*.tif*', 'Select TIF files', 'MultiSelect', 'off');
file = fullfile(path,filename);

% Read file, convert to double rep. and adjust contrast
image = im2double(imread(file));
image = imadjust(image);
%image = imcomplement(image);

% Plot image
figure(1)
imshow(image);

%% Parameters

px2mm = 0.08989; %1/20; %mm/px

improc = 0; % 1=process image (ah, medfilt, sharp)

%% Proccess image (if set)

if improc == 1
   
   %image = adapthisteq(image);
   image = medfilt2(image);     % Median filter to reduce noise
   image = imsharpen(image);    % Sharpening of the edges
    
end

%% Contrast-to-Noise Ratio (CNR)

image = imcomplement(image);    % Negative of image

%Define ROI for signal
sROI = drawrectangle();
sM = createMask(sROI);

avS = mean(image(sM));
stdS = std(image(sM));

%Define ROI for background
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
eROI = drawrectangle();
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
ylabel('')
legend('edge','ESF','LSF')

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
MTF = MTF/max(MTF);

%fX = linspace(0,length(xData)/2,numel(MTF)/2)/px2mm;

dx = (xData(2) - xData(1))*px2mm;   % Spatial step
df = 1/dx;                          % Spatial freq. step

fX = df*(0:length(xData)-1)/length(xData);     % Spatial freq. range
fX(fX>df/2) = fX(fX>df/2)-df;

figure(3)
plot(fX(1:15),MTF(1:15),'LineWidth', 2)
xlabel('cy/mm');
ylabel('MTF');

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