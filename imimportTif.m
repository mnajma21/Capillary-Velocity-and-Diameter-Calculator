function imageLines = imimportTif(fname)
% function imageLines = imimportTif(fname);
%
% This function imports a multiFrame TIF file as a single
% array. Modify this section to accomodate your raw data format.
%
% This particular file format assumes the time axis is along the
%   x axis (e.g. each y scan is instantaneous).  The output format
%   is spatial along the first axis, and time along the second axis.
%
%  Arbitrary limit of 256 frames in the TIF (which can be increased)
%
% PWG 3/28/2011


imageFrames = [];
disp('Frame: ');
for index = 1:512
    try
        imageFrames(:,:,index) = imread(fname,index);
    catch
        break
    end
    
    % give the user status updates
    fprintf('%d, ', index);
    if mod(index,10) == 0
        fprintf('\n')
    end
    if index == 512
        disp('Arbitrary 256 frame TIF limit reached. Increase if necessary')
    end
end
fprintf('\n');

% reformat into a single frame
imageFramesSize = size(imageFrames);

imageLines = zeros(imageFramesSize(2), imageFramesSize(1)*imageFramesSize(3));
for index = 1:imageFramesSize(3)
    imageLines(:,1+(index-1)*imageFramesSize(1):(index)*imageFramesSize(1)) = imageFrames(:,:,index)';
end

disp('Multiframe import complete');
