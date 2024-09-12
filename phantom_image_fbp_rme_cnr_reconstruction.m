% Save this as 'run_fbp_reconstruction.m'

% Clear environment
clc; clear; close all;

% Call the main function
fbp_reconstruction();

% Main Function
function fbp_reconstruction()
    % Define the size of the image
    image_size = 256; % Size of the square image to be reconstructed
    theta = 0:179;    % Projection angles (0 to 179 degrees)

    % Create a test object (Phantom)
    phantom_image = phantom(image_size);

    % Generate the sinogram (projections)
    projections = radon(phantom_image, theta);

    % Perform the filtered back projection with different filters
    filters = {'Ram-Lak', 'Shepp-Logan', 'Cosine', 'Hann', 'Flattop', 'Parzen'};
    
    % Preallocate arrays to store RMSE and CNR values for each filter
    rmse_values = zeros(1, length(filters));
    cnr_values = zeros(1, length(filters));
    
    figure;
    for i = 1:length(filters)
        filter_type = filters{i};
        % Filter and reconstruct the image
        reconstructed_image = filtered_back_projection(projections, theta, image_size, filter_type);
        
        % Calculate RMSE
        rmse_values(i) = calculate_rmse(reconstructed_image, phantom_image);
        
        % Calculate CNR
        cnr_values(i) = calculate_cnr(reconstructed_image);
        
        % Display the reconstructed image
        subplot(2, 3, i);
        imshow(reconstructed_image, []);
        title([filter_type, ', RMSE: ', num2str(rmse_values(i)), ', CNR: ', num2str(cnr_values(i))]);
    end
    
    % Display RMSE and CNR results
    disp('RMSE Values:');
    disp(array2table(rmse_values, 'VariableNames', filters));
    disp('CNR Values:');
    disp(array2table(cnr_values, 'VariableNames', filters));
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
            filter = freq;
            
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

function rmse = calculate_rmse(reconstructed_image, original_image)
    % Calculate the Root Mean Squared Error (RMSE)
    error = reconstructed_image - original_image;
    rmse = sqrt(mean(error(:).^2));
end

function cnr = calculate_cnr(reconstructed_image)
    % Calculate the Contrast-to-Noise Ratio (CNR)
    % Define object and background regions manually (this is an example)
    object_region = reconstructed_image(90:120, 90:120); % Example ROI for object
    background_region = reconstructed_image(200:230, 200:230); % Example ROI for background
    
    % Calculate mean and standard deviation for object and background
    mean_object = mean(object_region(:));
    std_object = std(object_region(:));
    mean_background = mean(background_region(:));
    std_background = std(background_region(:));
    
    % Compute CNR (Contrast-to-Noise Ratio)
    cnr = abs(mean_object - mean_background) / sqrt(std_object^2 + std_background^2);
end
