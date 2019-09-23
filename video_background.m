function [background, maxV, minV] = video_background(V, path_name, fn, frames, pixels_range, approach, ext, complement)

% Computes mean value of the whole movie to later on substract it from each
% frame. Substracting the background will get rid of all noise as well as
% of all particles that are not moving. It also calculates the limits of
% the de-backgrounded images (minV and maxV) for later histogram
% equalization
%
% input variables:
%
%       V is the videoreader or imfinfo object.
%
%       frames is a vector with the list of frames to be used for
%           thresolding.
%
%       pixels_range is a 2X2 matrix with matrix indexes for the
%          portion of the video frames to be used in the analyses.
%
%       Approach is the averaging method: 'mode', 'mean' or 'median'.
%           Default is median.
%
% output variables:
%
%       background is the average frame that is substracted from each frame
%       to get rid of non moving particles and noise
%
%       maxV and minV are the maximum and minimum intensity values of the
%           video after substracting the background and adding the mean
%           limits of the de-backgrounded images for histogram equalization
%
%   ex.:  [background, maxV, minV] = video_background(V, frames)
%
% Oscar Guadayol, 2015
% oscar_at_guadayol.cat
%
% Copyright (C)
%
% This file is part of trackbac.
%
% trackbac is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with trackbac.  If not, see <http://www.gnu.org/licenses/>.
%

% if strcmpi(ext,'.avi')
%     V = VideoReader(fullfile(path_name, [fn, ext]));
%     l = V.Duration*V.FrameRate;
% elseif strcmpi(ext, '.tif') ||strcmpi(ext, '.tiff')
%     V = imfinfo(fullfile(path_name, [fn, ext]));
% end

if isempty(approach)
    approach = 'median';
end

height = diff(pixels_range(1,:))+1;
width = diff(pixels_range(2,:))+1;

%% mode approach
% This works best when there are no pixels out of range (i.e., no strings
% or 0s or 255s
if strcmp(approach,'mode')
    if strcmpi(ext,'.avi')
        I = read(V);
        if strcmp(get(V,'videoFormat'),'RGB24')
            I = I(:,:,1,:)*0.2989+I(:,:,2,:)*0.5870+I(:,:,3,:)*0.1140;
            I = squeeze(I);
        end
    elseif strcmpi(ext,'.tif') || strcmpi(ext,'.tiff')
        I = zeros(V(1).Height, V(1).Width, length(frames), 'uint8');
        TifLink = Tiff(V(1).Filename, 'r');
        for ff=frames
            TifLink.setDirectory(ff);
            I(:,:,ff)=TifLink.read();
        end
        TifLink.close();
    elseif strcmpi(ext,'.mat')
        fns = dir([path_name fn(1:end-5),'*.png']);
        I = zeros(V(1).Height, V(1).Width, length(frames), 'uint8');
        parfor ff=1:length(frames)
            I(:,:,ff) = imread([path_name fns(ff).name]);
        end
    end
    I = I(pixels_range(1,1):pixels_range(1,2),pixels_range(2,1):pixels_range(2,2),frames);
    if complement
        I = imcomplement(I);
    end
    % background = mode(I,3); % this is very consuming
    background = zeros(height,width, 'uint8');
    parfor ii = 1:size(I,1)
        background(ii,:) = mode(I(ii,:,:),3);
    end
    maxV = max(I,[],3);
    minV = min(I,[],3);
    
    %% mean approach
elseif strcmp(approach,'mean')
    background = zeros(height,width');
    maxV = zeros(height,width);
    minV = ones(height,width)*255;
    
    if strcmpi(ext,'.avi')
        for ff = frames
            I = read(V,ff);
            if strcmp(get(V,'videoFormat'),'RGB24') % if there are 3 channels compute mean
                I = rgb2gray(I);
            end
            I = I(pixels_range(1,1):pixels_range(1,2),pixels_range(2,1):pixels_range(2,2));
            if complement
                I = imcomplement(I);
            end
            background = background + (1/length(frames)).*double(I);
            maxV = max(maxV, double(I));
            minV = min(minV, double(I));
            % background = background + (1/l).*mean(double(read(V,ii)),3);
        end
    elseif strcmpi(ext,'.tif')  || strcmpi(ext,'.tiff')
        for ff = frames
            I = imread(V(1).Filename,ff);
            I = I(pixels_range(1,1):pixels_range(1,2),pixels_range(2,1):pixels_range(2,2));
            if V(1).BitDepth>8
                if isempty(regexp(V(1).ImageDescription,'Hamamatsu', 'once')) && V(1).UnknownTags.Value==12
                    I = uint8(bitshift(I,-4)); % in fact, when the V.Bitdepth =16 in a hamamatsu multitif,the bitdepth=12
                else
                    I = uint8(I);
                end
            end
            if complement
                I = imcomplement(I);
            end
            background = background + (1/length(frames)).*double(I);
            maxV = max(maxV, double(I));
            minV = min(minV, double(I));
        end
    elseif strcmpi(ext,'.mat')
        fns = dir([path_name fn(1:end-5),'*.png']);
        for ff = frames
            I = imread([path_name fns(ff).name]);
            I = I(pixels_range(1,1):pixels_range(1,2),pixels_range(2,1):pixels_range(2,2));
            if complement
                I = imcomplement(I);
            end
            background = background + (1/length(frames)).*double(I);
            maxV = max(maxV, double(I));
            minV = min(minV, double(I));
        end
    end
    
    %% median approach
elseif strcmp(approach,'median')
    maxV = zeros(height,width);
    minV = ones(height,width)*255;
    if strcmpi(ext,'.avi')
        
        l = V.Duration*V.FrameRate;
        
        %% incremental median. It uses much less memory.
        % preallocates background
        I = zeros(V.Height, V.Width, 256, 'int16');
        V = VideoReader(fullfile(path_name, [fn, ext]));
        for ff = 1:l            
            frame = read(V, ff);
            if strcmp(get(V,'videoFormat'),'RGB24')
                frame = frame(:,:,1)*0.2989+frame(:,:,2)*0.5870+frame(:,:,3)*0.1140;
                frame = squeeze(frame);
            end
            if complement
                frame = imcomplement(frame);
            end
            
            ii = size(frame,1)*size(frame,2)*(double(frame(:)))+(1:numel(frame))';
            I(ii) = I(ii)+1;
            maxV = max(maxV, double(frame));
            minV = min(minV, double(frame));            
        end
        
    elseif  strcmpi(ext,'.tif')  || strcmpi(ext,'.tiff')
        if V(1).BitDepth==16
            bits = 'uint16';
        elseif V(1).BitDepth==8
            bits = 'uint8';
        end
        l = length(frames);
        % preallocates background
        I = zeros(V(1).Height, V(1).Width, 256, 'int16');
        TifLink = Tiff(V(1).Filename, 'r');
        for ff=frames
            TifLink.setDirectory(ff);
            frame = TifLink.read();
            if V(1).BitDepth==16
                if ~isempty(regexp(V(1).ImageDescription,'Hamamatsu', 'once')) && V(1).UnknownTags.Value==12
                    frame = uint8(bitshift(frame,-4)); % in fact, when the V.Bitdepth =16 in a hamamatsu multitif,the bitdepth=12
                else
                    frame = uint8(frame);
                end
            end
            if complement
                frame = imcomplement(frame);
            end
            ii = size(frame,1)*size(frame,2)*(double(frame(:)))+(1:numel(frame))';
            I(ii) = I(ii)+1;
            maxV = max(maxV, double(frame));
            minV = min(minV, double(frame));
        end
        TifLink.close();
        
    elseif  strcmpi(ext,'.mat')
        fns = dir([path_name fn(1:end-5),'*.png']);
        l = length(frames);
        I = zeros(height,width,length(frames), 'int8');
        for ff = frames
            frame = imread([path_name fns(ff).name]);
            ii = size(frame,1)*size(frame,2)*(double(frame(:)))+(1:numel(frame))';
            I(ii) = I(ii)+1;
            maxV = max(maxV, double(frame));
            minV = min(minV, double(frame));
        end
        if complement
            frame = imcomplement(frame);
        end
    end
    
    t = cumsum(I,3);
    clear I
    t = t>=l/2;
    t = uint8(t);
    t = diff(t,1,3);
    [r,c,k] = ind2sub(size(t),find(t));
    clear t
    background = uint8(zeros(size(frame)));
    background(sub2ind(size(frame),r,c)) = k;    
end

maxV = max(max((double(maxV)-double(background)+min(double(background(:))))));
minV = min(min((double(minV)-double(background)+min(double(background(:))))));
background = uint8(background);
