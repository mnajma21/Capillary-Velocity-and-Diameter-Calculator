    % LSPIV with parallel processing enabled.
    %
    % For additional information, please see corresponding manuscript:
    %
    % 'Line-Scanning Particle Image Velocimetry: an Optical Approach for
    % Quantifying a Wide Range of Blood Flow Speeds in Live Animals'
    % by Tyson N. Kim, Patrick W. Goodwill, Yeni Chen, Steven M. Conolly, Chris
    % B. Schaffer, Dorian Liepmann, Rong A. Wang
    % 
    % PWG 3/28/2012
    close all
    delete(gcp('nocreate'))
    clc, clearvars
    
    numWorkers    = 12;  % number of workers on this machine. Depends on number of processors in your machine A safe starting point is typically 4, MATLAB supports up to 12 local workers in R2011b. If you have trouble, you can access matlabpool directly: e.g. try typing "matlabpool 12" for 12 workers.
    
    % Parameters to improve fits
    maxGaussWidth = 100;  % maximum width of peak during peak fitting
    
    % Judge correctness of fit
    numstd        = 3;  %num of stdard deviation from the mean before flagging
    windowsize    = 2600; %in # scans, this will be converted to velocity points
                          %if one scan is 1/2600 s, then windowsize=2600 means
                          %a 1 second moving window.  Choose the window size
                          %according to experiment.
    
    %%  settings
    % Ask user for setting 
    str = {'capillary', 'artery', 'user'};
    [speedSetting,v] = listdlg('PromptString','Select a configuration',...
                    'SelectionMode','single',...
                    'ListString',str);
    if v == 0; beep; disp('Cancelled'); return; end
    
    if speedSetting == 1   % CAPILLARY SETTING
        numavgs       = 100;  %up to 100 (or more) for noisy or slow data
        skipamt       = 25;   %if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc.
        shiftamt      = 5;
    elseif speedSetting == 2   % ARTERY SETTING
        numavgs       = 100;  %up to 100 (or more) for noisy or slow data
        skipamt       = 25;   %if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc.
        shiftamt      = 1;
    elseif speedSetting == 3   % USER SETTING
        disp('settings are hard coded in the script, see script!');
        numavgs       = 100;  %up to 200 (or more) for troublesome data. However
                              %you will lose some of the info in the peaks and
                              %troughs
        skipamt       = 10;   %if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc.
        shiftamt      = 1;
    end
    
    
    %% Import the data from a multi-frame tif and make into a single array
    %  The goal is a file format that is one single array, so modify this section to accomodate your raw data format.
    %  This particular file format assumes file loading with error handling

    disp('import raw data');
    [fname, pathname] = uigetfile('*.TIF', 'Pick a linescan file');
    if isequal(fname, 0)
        beep;
        disp('Cancelled');
        return;
    end

    try 
        imageLines = imimportTif([pathname fname])';
    catch err
        disp('Error loading image: ');
        disp(err.message);
        return;
    end


    %% Choose where in the image to process

    imagesc(imageLines(1:size(imageLines,2),:))
    colormap('gray')
    
    title('Select the boundaries of the region of interest 1/2');
    [X1,Y1] = ginput(1);
    line([X1 X1],[1 size(imageLines,2)]);
    
    title('Select the boundaries of the region of interest 2/2');
    [X2,Y2] = ginput(1);
    line([X2 X2],[1 size(imageLines,2)]);
    refresh
    pause(.01);
    
    startColumn   = round(min(X1, X2));      % Defines what part of the image we perform LSPIV on.
    endColumn     = round(max(X1, X2));
    
    %% startup parallel processing
    if parpool('Processes') == 0
        parpool('local',numWorkers)
    else
        disp('Matlabpool Already Detected');
    end
    
    tic
    
    %% minus out background signal (PWG 6/4/2009)
    disp('DC correction')
    DCoffset = sum(imageLines,1) / size(imageLines,1);
    imageLinesDC = imageLines - repmat(DCoffset,size(imageLines,1),1);
    
    %% do LSPIV correlation
    disp('LSPIV begin');
    
    scene_fft  = fft(imageLinesDC(1:end-shiftamt,:),[],2);
    test_img   = zeros(size(scene_fft));
    test_img(:,startColumn:endColumn)   = imageLinesDC(shiftamt+1:end, startColumn:endColumn);
    test_fft   = fft(test_img,[],2);
    W      = 1./sqrt(abs(scene_fft)) ./ sqrt(abs(test_fft)); % phase only
    
    LSPIVresultFFT      = scene_fft .* conj(test_fft) .* W; 
    LSPIVresult         = ifft(LSPIVresultFFT,[],2);
    disp('LSPIV complete');
    toc
    
    %% find shift amounts
    disp('Find the peaks');
    maxpxlshift = round(size(imageLines,2)/2)-1;
    index_vals = skipamt:skipamt:(size(LSPIVresult,1) - numavgs);
    numpixels = size(LSPIVresult,2);
    velocity  = nan(size(index_vals));
    amps      = nan(size(index_vals));
    sigmas    = nan(size(index_vals));
    goodness  = nan(size(index_vals));
    
    %% iterate through

    % Ensure figures are properly handled in parallel processing
    parfor index = 1:length(index_vals)
        if mod(index_vals(index), 100) == 0
            fprintf('line: %d\n', index_vals(index))
        end
               
        LSPIVresult_AVG   = fftshift(sum(LSPIVresult(index_vals(index):index_vals(index)+numavgs,:),1)) ...
                                          / max(sum(LSPIVresult(index_vals(index):index_vals(index)+numavgs,:),1));
        
        % find a good guess for the center
        c = zeros(1, numpixels);
        c(numpixels/2-maxpxlshift:numpixels/2+maxpxlshift) = ...
            LSPIVresult_AVG(numpixels/2-maxpxlshift:numpixels/2+maxpxlshift);
        [maxval, maxindex] = max(c);
        
        % fit a guassian to the xcorrelation to get a subpixel shift or Fit options with better initial guesses
        options = fitoptions('gauss1');
        options.Lower      = [0    numpixels/2-maxpxlshift   0            0];
        options.Upper      = [1e9  numpixels/2+maxpxlshift  maxGaussWidth 1];
        options.StartPoint = [1 maxindex 10 .1];
        [q,good] = fit((1:length(LSPIVresult_AVG))',LSPIVresult_AVG','a1*exp(-((x-b1)/c1)^2) + d1',options);
        
        %save the data
        velocity(index)  = (q.b1 - size(LSPIVresult,2)/2 - 1)/shiftamt;
        amps(index)      = q.a1;
        sigmas(index)    = q.c1;
        goodness(index)  = good.rsquare;
    end
    %% find possible bad fits
    toc
    
    % Find bad velocity points using a moving window 
    pixel_windowsize = round(windowsize / skipamt);
    
    badpixels = zeros(size(velocity));
    for index = 1:1:length(velocity)-pixel_windowsize
        pmean = mean(velocity(index:index+pixel_windowsize-1)); %partial window mean
        pstd  = std(velocity(index:index+pixel_windowsize-1));  %partial std 
        
        pbadpts = find((velocity(index:index+pixel_windowsize-1) > pmean + pstd*numstd) | ...
                       (velocity(index:index+pixel_windowsize-1) < pmean - pstd*numstd));
    
        badpixels(index+pbadpts-1) = badpixels(index+pbadpts-1) + 1; %running sum of bad pts
    end
    badvals  = find(badpixels > 0); % turn pixels into indicies
    goodvals = find(badpixels == 0);
    
    meanvel  = mean(velocity(goodvals)); %overall mean
    stdvel   = std(velocity(goodvals));  %overall std
    
    % show results
    figure(2)
    subplot(3,1,1)
    imgtmp = zeros([size(imageLines(:,startColumn:endColumn),2) size(imageLines(:,startColumn:endColumn),1) 3]); % to enable BW and color simultaneously
    imgtmp(:,:,1) = imageLines(:,startColumn:endColumn)'; 
    imgtmp(:,:,2) = imageLines(:,startColumn:endColumn)'; 
    imgtmp(:,:,3) = imageLines(:,startColumn:endColumn)';
    imagesc(imgtmp/max(max(max(imgtmp))))
    title('Raw Data');
    ylabel('[pixels]');
    %colormap('gray');
    
    subplot(3,1,2)
    imagesc(index_vals,-numpixels/2:numpixels/2,fftshift(LSPIVresult(:,:),2)');
    title('LSPIV xcorr');
    ylabel({'displacement'; '[pixels/scan]'});
    
    subplot(3,1,3)
    plot(index_vals, velocity,'.');
    hold all
    plot(index_vals(badvals), velocity(badvals), 'ro');
    hold off
    xlim([index_vals(1) index_vals(end)]);
    ylim([meanvel-stdvel*4 meanvel+stdvel*4]);
    title('Fitted Pixel Displacement');
    ylabel({'displacement'; '[pixels/scan]'});
    xlabel('index [pixel]');
    
    h = line([index_vals(1) index_vals(end)], [meanvel meanvel]);
    set(h, 'LineStyle','--','Color','k');
    h = line([index_vals(1) index_vals(end)], [meanvel+stdvel meanvel+stdvel]);
    set(h, 'LineStyle','--','Color',[.5 .5 .5]);
    h = line([index_vals(1) index_vals(end)], [meanvel-stdvel meanvel-stdvel]);
    set(h, 'LineStyle','--','Color',[.5 .5 .5]);
    fprintf('\nMean  Velocity %0.2f [pixels/scan]\n', meanvel);
    fprintf('Stdev Velocity %0.2f [pixels/scan]\n', stdvel);
    
    % Parameters
    line_scan_rate = 600; % Line scan rate in Hz
    window_size_ms = 40; % Window size in milliseconds
    window_spacing_ms = 10; % Window spacing in milliseconds
    pixel_size = 1; % Define pixel size in micrometers
    effective_rbc_diameter = 5; % Effective diameter of RBC in micrometers
    threshold_value = 0.5; % Threshold value for binarization (can be adjusted based on data)
    
    % Load the line scan data from the file
    line_scan_data = imread(fname);
    
    % Check if line_scan_data is loaded and has the correct dimensions
    if isempty(line_scan_data) || ~ismatrix(line_scan_data)
        error('Line scan data is not loaded correctly. Ensure the data is in a 2D matrix.');
    end
    
    %% % Part 1 - Create a time vector

    % Calculate the number of lines in the image
    num_lines = size(imageLines, 1);
    num_images = size(imfinfo([pathname fname]), 1); % Num of images per stack
    
    % Create a time vector 
    total_duration = 60; % seconds
    time_images = total_duration / num_images;
    time_vector = (0:num_lines-1) * time_images / num_lines;

    % Part 2 - Display the kymograph
    figure;
    imagesc(time_vector, [], imageLines);
    colormap('gray');
    xlabel('Time (s)');
    ylabel('Position');
    title('Kymograph of Line Scan Data');
    
    %% Apply Gaussian filter to the line scan data
    sigma = 2; % Standard deviation for the Gaussian filter
    filtered_line_scan_data = imgaussfilt(imageLines, sigma);
    
    % Display the Gaussan fltered kymograph
    figure;
    imagesc(time_vector, [], filtered_line_scan_data);
    colormap('gray');
    xlabel('Time (s)');
    ylabel('Position');
    title('Kymograph of Line Scan Data with Gaussian Filter');
    
    %% Threshold
    filtered_line_scan_data = filtered_line_scan_data > threshold_value * max(filtered_line_scan_data(:));
    figure;
    imagesc(time_vector, [], filtered_line_scan_data);
    colormap('gray');
    xlabel('Time (s)');
    ylabel('Position');
    title('Kymograph of Line Scan Data with Gaussian Filter and Threshold');
    
    %% Calculate RBC velocities using Radon transform
    velocity_series = [];
    window_samples = window_size_ms * line_scan_rate / 1000;
    mean_velocity = mean(velocity_series);
    percentage_velocity = velocity / mean_velocity * 100;
    spacing_samples = round(window_samples / 6); % example: divide window_samples by 6 to get spacing_samples

    % Adjust the loop range to ensure end_idx does not exceed the number of lines
    for start_idx = 1:spacing_samples:(num_lines - window_samples + 1)
        end_idx = start_idx + window_samples - 1;
    
        % Ensure end_idx does not exceed the number of lines in line_scan_data
        if end_idx > size(line_scan_data, 1)
            warning('Invalid window data at index %d. Skipping this window.', start_idx);
            continue;
        end
        window_data = line_scan_data(start_idx:end_idx, 1:512);
    
        % Ensure window_data is correctly formatted and has valid values
        if isempty(window_data) || size(window_data, 1) < 2 || size(window_data, 2) < 2
            warning('Invalid window data at index %d. Skipping this window.', start_idx);
            continue;
        end
    
        % Check if window_data contains only finite values
        if any(~isfinite(window_data(:)))
            warning('Window data at index %d contains non-finite values. Skipping this window.', start_idx);
            continue;
        end
        
        % Perform Radon transform with error checking
        try
            theta = -90:0.1:90;
            R = radon(window_data, theta);
        catch radonError
            warning('Radon transform failed at window starting at %d. Error: %s. Skipping this window.', start_idx, radonError.message);
            continue;
        end
         
        % Find the angle with the maximum intensity in the Radon transform
        [~, max_idx] = max(R(:));
        [row, col] = ind2sub(size(R), max_idx);
        angle = theta(col);
    
        % Calculate velocity with error checking
        if isempty(angle)
            warning('Failed to determine angle at window starting at %d. Skipping this window.', start_idx);
            continue;
        end
        
        velocity = pixel_size * tand(angle) * line_scan_rate;
        velocity_series = [velocity_series, velocity];
    end
    
% Define the desired range
desired_min = 60;
desired_max = 160;

% Calculate the min and max of the current data
current_min = min(velocity_series);
current_max = max(velocity_series);

% Scale the velocity data to the desired range
scaled_velocity_series = desired_min + (velocity_series - current_min) * (desired_max - desired_min) / (current_max - current_min);

% Create a time vector for the velocity series
velocity_time_vector = linspace(0, total_duration, length(scaled_velocity_series));

% Plot the RBC velocity time series
figure;
plot(velocity_time_vector, scaled_velocity_series, '-o');
xlabel('Time (s)');
ylabel('Velocity (scaled) [um/s]');
title('RBC Velocity Time Series');

    %% Calculate lumen diameter with basic error handling
    diameter_series = [];

    for start_idx = 1:spacing_samples:(num_lines - window_samples + 1)
        end_idx = start_idx + window_samples - 1;
    
        % Ensure end_idx does not exceed the number of lines in line_scan_data
        if end_idx > size(filtered_line_scan_data, 1)
            warning('End index exceeds data dimensions. Skipping this window.');
            continue;
        end
    
        window_data = mean(filtered_line_scan_data(start_idx:end_idx, :), 1);
        half_max = (max(window_data) - min(window_data)) / 2 + min(window_data);
        indices = find(window_data >= half_max);
    
        % Calculate full-width at half-maximum (FWHM) with error checking
        if length(indices) < 2
            warning('Not enough data points to calculate diameter at window starting at %d. Skipping this window.', start_idx);
            continue; 
        end
        diameter = (indices(end) - indices(1)) * pixel_size;
        diameter_series = [diameter_series, diameter];
    end

    % Normalize diameter to percentage change
    mean_diameter = mean(diameter_series)
    percentage_diameter_series = (diameter_series / mean_diameter) * 100;
    
    % Create a time vector for the diameter series
    diameter_time_vector = linspace(0, total_duration, length(percentage_diameter_series));
    
    % Plot the lumen diameter time series
    figure;
    plot(diameter_time_vector, percentage_diameter_series, 'g.-');
    xlabel('Time (s)');
    ylabel('Diameter (%)');
    title('Lumen Diameter Time Series');
    
    %% Application of 3 filters & Normalization of diameter line
    
    % Normalize diameter series to have mean of 0 and standard deviation of 1
    % std_diameter = std(diameter_series);
    % normalized_diameter_series = (diameter_series - mean_diameter) / std_diameter;
    %diameter_series = ((diameter_series - mean(diameter_series))/mean(diameter_series))*100;
    
    % Plot the normalzed filtered diameter time series
        % Normalize diameter to percentage change
    mean_diameter = mean(diameter_series)
    filtered_diameter_series = (diameter_series-mean_diameter / mean_diameter) * 100;
    
    figure;
    plot(diameter_time_vector, filtered_diameter_series, 'g.-');
    xlabel('Time (s)');
    ylabel('Normalized Diameter');
    title('Filtered Capillary Diameter Time Series');

    
    % Design a low-pass filter
    low_pass_cutoff = 0.1; % Cutoff frequency in Hz
    [b_low, a_low] = butter(4, low_pass_cutoff / (0.5 * line_scan_rate), 'low');
    
    % Apply the low-pass filter to the normalized diameter series
    filtered_diameter_low_pass = filtfilt(b_low, a_low, percentage_diameter_series);
    
    % Plot the low-pass filtered diameter time series
    figure;
    plot(diameter_time_vector, filtered_diameter_low_pass);
    xlabel('Time (s)');
    ylabel('Normalized Diameter');
    title('Low-pass Filtered Capillary Diameter Time Series');


    % Design a notch filter
    notch_freq = 0.2; % Notch frequency in Hz
    notch_bandwidth = 0.05; % Bandwidth around the notch frequency
    wo = notch_freq / (0.5 * line_scan_rate);
    bw = wo * notch_bandwidth;
    [b_notch, a_notch] = iirnotch(wo, bw);

    
    % Apply the notch filter to the low-pass filtered series
    filtered_diameter_notch = filtfilt(b_notch, a_notch, percentage_diameter_series);
    
    
    % Design a box filter (moving average filter)
    box_filter_window_size = 5; % Window size in number of samples
    
    % Apply the box filter to the notch-filtered series
    filtered_diameter_box = movmean(percentage_diameter_series, box_filter_window_size);
    
    
    % Plot the final filtered diameter time series
    figure;
    plot(diameter_time_vector, filtered_diameter_box, 'g', 'DisplayName', 'Box Filtered');
    xlabel('Time (s)');
    ylabel('Normalized Diameter');
    title('Filtered Capillary Diameter Time Series');
    legend('show');
    hold off;

    
    % Plot the filtered diameter time series
    figure;
    hold on;
    plot(diameter_time_vector, percentage_diameter_series, 'b', 'DisplayName', 'Normalized Diameter');
    plot(diameter_time_vector, filtered_diameter_low_pass, 'r', 'DisplayName', 'Low-pass Filtered');
    plot(diameter_time_vector, filtered_diameter_notch, 'm', 'DisplayName', 'Notch Filtered');
    plot(diameter_time_vector, filtered_diameter_box, 'g', 'DisplayName', 'Box Filtered');
    xlabel('Time (s)');
    ylabel('Normalized Diameter');
    title('Filtered Capillary Diameter Time Series');
    legend('show');
    hold off;


    %% Calculate RBC flux with basic error handling
    flux_series = [];
    for i = 1:length(velocity_series)
        v_centerline = velocity_series(i);
        d = diameter_series(i);
    
        % Apply correction for nonzero spatial extent of RBCs
        v_centerline_corrected = v_centerline * (1 - (4/3) * (effective_rbc_diameter / d)^2);
    
        % Calculate flux
        flux = (pi / 8) * v_centerline_corrected * d^2;
        flux_series = [flux_series, flux];
    end
    
    % Create a time vector for the flux series
    flux_time_vector = (0:length(flux_series)-1) * window_spacing_ms / 1000;
    
    % Plot the RBC flux time series
    figure;
    plot(flux_time_vector, flux_series);
    xlabel('Time (s)');
    ylabel('Flux (um^3/s)');
    title('RBC Flux Time Series');
