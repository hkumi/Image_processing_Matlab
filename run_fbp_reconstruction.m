
% Clear environment
clc; clear; close all;

% Main Function
fbp_reconstruction();

% Clear environment
clc; clear; close all;

% Main Function
fbp_reconstruction();

function fbp_reconstruction()
    % Define the size of the image
    image_size = 256; % Size of the square image to be reconstructed
    theta = 0:179;    % Projection angles (0 to 179 degrees)

    % Load a custom image (change this to your image file path)
    addpath(genpath(pwd))
    custom_image = imread('Cone-Beam-and-3D-Imaging.jpg.optimal.jpg'); % Replace with your image file
    custom_image = imresize(custom_image, [image_size image_size]); % Resize to fit
    custom_image = rgb2gray(custom_image); % Convert to grayscale if it's RGB

    % Generate the sinogram (projections) from the custom image
    projections = radon(custom_image, theta);

    % Perform the filtered back projection with different filters
    filters = {'Ram-Lak', 'Shepp-Logan', 'Cosine', 'Hann', 'Flattop', 'Parzen'};
    
    figure;
    for i = 1:length(filters)
        filter_type = filters{i};
        % Filter and reconstruct the image
        reconstructed_image = filtered_back_projection(projections, theta, image_size, filter_type);
        
        % Display the reconstructed image
        subplot(2, 3, i);
        imshow(reconstructed_image, []);
        title(filter_type);
    end
end
function reconstructed_image = filtered_back_projection(projections, theta, image_size, filter_type)
    % Get the size of projections
    [num_projections, num_angles] = size(projections);
    
    % Frequency range for the filter
    filter = generate_filter(filter_type, num_projections);

    % Apply the filter in the frequency domain
    for i = 1:num_angles
        projection_fft = fft(projections(:, i));  % FFT of the projection
        % Zero-padding the filter to match the length of projection_fft
        padded_filter = ifftshift(filter); % Apply shift to center
        projection_fft = projection_fft .* padded_filter; % Apply the filter
        projections(:, i) = real(ifft(projection_fft)); % Inverse FFT
    end

    % Perform the back projection
    reconstructed_image = iradon(projections, theta, 'linear', 'none', 1, image_size);
end

function filter = generate_filter(filter_type, num_projections)
    % Generate the filter based on the input type
    n = ceil(num_projections/2); % Half length for symmetry
    freq = linspace(0, 1, n)'; % Frequency range from 0 to 1
    
    switch filter_type
        case 'Ram-Lak'
            % Ram-Lak (ramp) filter
            filter = abs(freq);
            
        case 'Shepp-Logan'
            % Shepp-Logan filter
            filter = freq .* (sin(pi * freq) ./ (pi * freq));
            filter(1) = 1; % Handle division by zero at freq = 0
            
        case 'Cosine'
            % Cosine filter
            filter = freq .* cos(pi * freq / 2);
            
        case 'Hann'
            % Hann filter
            filter = freq .* (0.5 + 0.5 * cos(pi * freq));
            
        case 'Flattop'
            % Flattop filter
            a0 = 0.21557895; a1 = 0.41663158; a2 = 0.277263158;
            a3 = 0.083578947; a4 = 0.006947368;
            filter = a0 - a1 * cos(2 * pi * freq) + ...
                     a2 * cos(4 * pi * freq) - ...
                     a3 * cos(6 * pi * freq) + ...
                     a4 * cos(8 * pi * freq);
                 
        case 'Parzen'
            % Parzen filter (triangular filter)
            filter = freq .* (1 - abs(freq / max(freq)));
            
        otherwise
            % Default to Ram-Lak (no additional modifications)
            filter = freq;
    end

    % Symmetrize the filter for both positive and negative frequencies
    filter = [filter; flipud(filter(1:end-1))];
end
