function [Boundaries, Centroids, Geometries, fn, boundary_threshold, ...
    particle_threshold, frames, background, minV, maxV, FrameRate] = ...
    particle_detection(fn, path_name, varargin)

% This function automatically detects moving particles in a video file (avi
% of stacked tiff) and records geometrical parameters of each particle at
% each particular frame. It then saves the results in a '.mat' file
% with the same name and location as the original video.
%
% input variables:
%
%       fn is the filename (with or without path, optional). If empty, it
%           opens a dialog to select file.
%
%       path_name is the pathname (optional). If empty, and fn has no path,
%           it opens a dialog to select file.
%
%       scaling in pixels/microm (optional). If no input , data will be
%           expressed in pixels;
%
%       boundary_threshold (optional) is the minimum pixel intensity in a
%           grayscale (0-255) of the boundaries of a particle. If no value
%           is given, the function threshold_level function is called.
%
%       particle_threshold (optional) is the minimum pixel intensity in a
%           grayscale (0-255) that a particle needs to have to be included
%           in the analyses. If no value is given, the function
%           threshold_level function is called.
%
%       frame_range (optional) is the first and last frames to be used by
%           the script. This is useful when the avi file is corrupted at
%           the beggining or the end of the file. Default is the complete
%           video.
%
%       pixels_range (optional) is a 2X2 matrix with matrix indexes for the
%          portion of the video frames to be used in the analyses. Useful
%          to cut off information.
%
%       FrameRate (optional) is the frame rate of the video. If empty,
%           information is retrieved from the avi.
%
%       background is the average frame that is substracted from each frame
%           to get rid of non moving particles and noise.
%
%       maxV and minV are the maximum and minimum intensity values of the
%           video after substracting the background and adding the mean
%           limits of the de-backgrounded images for histogram
%           equalization.
%
%       complement is a logical variable. If true, images are inverted
%           using imcomplement. Useful for light particles on dark field,
%           such as fluorescent images. Default is false.
%
%       noholes is a logic variable. If true, bwboundaries searches for
%           object (parent and child) boundaries. This can provide better
%           performance. If false, bwboundaries searches for both object
%           and hole boundaries. Default is false.
%
%       minimum_radius is the radius in µm of the smallest particle to be
%           considered a cell. Default is 0.5 µm
%
% output variables:
%
%       Boundaries is a cell array of size (NumberofFramesX1) in which each
%           cell is another array of 2-column matrices, one for each
%           particle in the frame, with x,y points of
%           the boundary of the particle.
%
%       Centroids is a cell array of size (NumberofFramesX1), in which each
%          cell is a NX2 matrix, N being the number of particles
%          detected within the frame, and the 2 columns represent position
%          x and y in pixelds of the centroids of the particles.
%
%       Geometries is a cell array of size (NumberofFramesX1), in which
%          each cell is a NX11 matrix, N being the number of particles.
%          Columns are area, MajorAxisLength, MinorAxisLength,
%          eccentricity, equivDiameter, orientation, perimeter and solidity
%          of each particle as defined in function regionprops, and
%          arc_length, minimum distance between poles, and width. Units are
%          pix
%
%       fn is the filename.
%
%       boundary_threshold (optional) is the minimum pixel intensity in a
%           grayscale (0-255) of the boundaries of a particle.
%
%       particle_thresholds (optional) is the minimum pixel intensity in
%           a grayscale (0-255) that a particle needs to have to be
%           included in the analyses.
%
%       frames is a vector with the list of frames to be used for
%           thresolding.
%
%       background is the average frame that is substracted from each frame
%           to get rid of non moving particles and noise.
%
%       maxV and minV are the maximum and minimum intensity values of the
%           video after substracting the background and adding the mean
%           limits of the de-backgrounded images for histogram
%           equalization.
%
%       FrameRate is the video frame rate in frames/s
%
%
% Requires: image processing toolbox, parfor_progress
%
% Calls: video_background, threshold_level, preprocess,
%   bwboundaries_noholes, Progress
%
% ex.:
%    [Boundaries, Centroids, Geometries, fn, boundary_threshold,...
%     particle_threshold, frames, background, minV, maxV, FrameRate] =...
%     particle_detection(fn, '', 'scaling',...
%       2.75, 'boundary_threshold', 100, 'frame_range', [1,100])
%
% 
%  Copyright (C) 2016,  Oscar Guadayol
%  oscar_at_guadayol.cat
%
%
%  This program is free software; you can redistribute it and/or modify it
%  under the terms of the GNU General Public License, version 3.0, as
%  published by the Free Software Foundation.
% 
%  This program is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%  General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License along
%  with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%  This files is part of trackbac, a collection of matlab scripts to
%  geometrically characterize and track swimming bacteria imaged by a
%  phase-contrast microscope

%% warnings
warning('error', 'MATLAB:rankDeficientMatrix');
warning('error', 'MATLAB:singularMatrix')
warning('off','MATLAB:inpolygon:ModelingWorldLower')
spmd % apply warning supressions to all workers
    warning('error', 'MATLAB:rankDeficientMatrix');
    warning('error', 'MATLAB:singularMatrix')
    warning('off','MATLAB:inpolygon:ModelingWorldLower')
end

%% Parses input arguments and loads data %%

p = inputParser;
addParameter(p,'scaling',[],@isnumeric);
addParameter(p,'boundary_threshold',[]);
addParameter(p,'particle_threshold',[]);
addParameter(p,'frame_range',[],@isnumeric);
addParameter(p,'pixels_range',[],@isnumeric);
addParameter(p,'FrameRate',[],@isnumeric);
addParameter(p,'background',[],@isnumeric);
addParameter(p,'maxV',[],@isnumeric);
addParameter(p,'minV',[],@isnumeric);
addParameter(p,'complement',[],@islogical);
addParameter(p,'noholes',[],@islogical);
addParameter(p,'minimum_radius',[],@isnumeric);

parse(p,varargin{:});
scaling = p.Results.scaling;

if isempty(p.Results.boundary_threshold)
    boundary_threshold = [];
elseif numel(p.Results.boundary_threshold)>0
    boundary_threshold = p.Results.boundary_threshold(1);
end

if isempty(p.Results.particle_threshold)
    particle_threshold = boundary_threshold;
elseif numel(p.Results.particle_threshold)>0
    particle_threshold = p.Results.particle_threshold;
end

FrameRate = p.Results.FrameRate;
background = p.Results.background;
minV = p.Results.minV;
maxV = p.Results.maxV;

if isempty(p.Results.complement)
    complement = false;
else
    complement = p.Results.complement;
end

if isempty(p.Results.noholes)
    noholes = false;
else
    noholes = p.Results.noholes;
end

if isempty(p.Results.minimum_radius)
    minimum_radius = 0.5;
else
    minimum_radius = p.Results.minimum_radius;
end

% selects video file
if isempty(fn)
    if isempty(path_name)
        path_name = pwd;
    end
    [fn, path_name] = uigetfile('*.avi; *.tif; *.tif; *info.mat', 'Choose an avi, tiff or mat file', path_name);
    fn = strcat(path_name,fn);    
    [~,fn] = fileattrib(fn);
    fn = fn.Name;
end

% changes \ to / to make it readable by fileparts in all systems
fn = strrep(fn, '\', filesep);
 
if ~isempty(path_name)
    [~, fn, ext] = fileparts(fn);
else
    [path_name, fn, ext] = fileparts(fn);
end
if strcmp(ext,'.mat')
    outfile = fullfile(path_name,[fn(1:end-5) '.mat']);
else
    outfile = fullfile(path_name,[fn '.mat']);
end
if isempty(scaling)
    warning('No scaling input. Data will be expressed in pixels')
    scaling = 1;
end

if strcmpi(ext,'.avi')
    V = VideoReader(fullfile(path_name, [fn, ext]))
    l = V.Duration*V.FrameRate;
    if isempty(FrameRate)
        FrameRate = get(V,'FrameRate');
    end
    StartTime = V.CurrentTime;
    if V.CurrentTime==0
        StartTime = dir(fullfile(path_name, [fn, ext]));
        StartTime = StartTime.datenum;
    end
elseif strcmpi(ext, '.tif') ||strcmpi(ext, '.tiff') 
    V = imfinfo(fullfile(path_name, [fn, ext]))
    
    % retrieves time information from multitiff captured using HCImage Live
    % software. Other software may need a different aproach to get the
    % times and frame rates.
 
    if ~isempty(regexp(V(1).ImageDescription,'Hamamatsu', 'once'))
        ImageDescription = char(V.ImageDescription);
        StartTime = datenum(V(1).ImageDescription(33:52));
        [~,tt0] = regexp(ImageDescription(1,:),'Time_From_Start = ');
        time = ImageDescription(:,tt0+1:tt0+13);
        time = datenum(time);
        if isempty(FrameRate)
            FrameRate = 1/mean(diff(time))/24/3600;
        end
    else
        StartTime = datenum(V(1).FileModDate);
        if isempty(FrameRate)
            FrameRate = 1;
        end
        time = (0:size(V,1)-1)/FrameRate;
    end
    l = size(V,1); % number of frames in video
    
elseif strcmpi(ext, '.mat')
    path_name = [path_name, filesep];
    load([path_name,fn])
    StartTime = datenum(V.ImageDescription(33:52));
    time = V.time;
    l = size(V.time,1);
    FrameRate = 1./mean(diff(time));
end

%% selects frames in the given range
if isempty(p.Results.frame_range)
    frames = 1:l;
else
    frames =  p.Results.frame_range(1):p.Results.frame_range(2);
end

%% selects pixels in each frame
if isempty(p.Results.pixels_range)
    pixels_range = [1,V(1).Height;1,V(1).Width];
else
    pixels_range =  p.Results.pixels_range;
end

%% minimum area, in pixels, a bacteria of 1µm radius would have
minimum_area = round((pi*(minimum_radius*scaling)^2))*2;

%% Computes mean value of the whole movie to later on substract it from each frame.
% This gets rid of all noise as well as of all particles that are not
% moving. It also calculates the limits of the de-backgrounded images for
% later histogram equalization
if any([isempty(background), isempty(maxV), isempty(minV)])
    [background, maxV, minV] = video_background(V, path_name, fn, frames, pixels_range, 'median',ext, complement);
end

save(outfile,'ext', 'fn', 'frames',...
    'background', 'minV',...
    'maxV', 'FrameRate','scaling', 'pixels_range', 'StartTime','-v7.3')

%% Determination of thresholds
if all([isempty(boundary_threshold),isempty(boundary_threshold)]) ||...
        all([isnan(boundary_threshold),isnan(boundary_threshold)])
    [boundary_threshold, particle_threshold] = threshold_level(V,...
        background, maxV, minV, frames, minimum_area, pixels_range,...
        scaling, complement, noholes, ext);
end
save(outfile,'boundary_threshold', 'particle_threshold', 'StartTime','-append')

%% Finds particles in each frame, and calculates geometrical parameters of each one
Boundaries = cell(l,1);
Centroids = cell(l,1);
Geometries = cell(l,1);

if ~isnan(boundary_threshold)
    Boundaries = cell(length(frames),1);
    Centroids = cell(length(frames),1);
    Geometries = cell(length(frames),1);
    N = length(frames);
    
    if isempty(gcp('nocreate'))
        c = parcluster;
        parpool(c.NumWorkers);
    end

    if strcmpi(ext,'.avi')
        parfor_progress(N);
        parfor ff = frames % frame
            V = VideoReader(fullfile(path_name, [fn, ext]));
            I = [];
            I = read(V,ff);
            if complement
                I = imcomplement(I);
            end
            I = I(pixels_range(1,1):pixels_range(1,2),...
                pixels_range(2,1):pixels_range(2,2));
            [I, ~] = image_normalization(I, background, minV, maxV, scaling);
            if ~isnan(boundary_threshold)
                % identifies only particles without holes
                [Boundaries{ff}, Centroids{ff}, Geometries{ff}]=...
                    bwboundaries_noholes(I, boundary_threshold,...
                    particle_threshold, minimum_area, false);
            end
            parfor_progress;
        end
        
    elseif strcmpi(ext, '.tif') ||strcmpi(ext, '.tiff') 
        parfor_progress(N);
        parfor ff = 1:frames(end) % frame
            I = [];
            I = imread(V(1).Filename,ff);
            if V(1).BitDepth==16
                if ~isempty(regexp(V(1).ImageDescription,'Hamamatsu', 'once')) || V(1).BitDepth==12
                    I = uint8(bitshift(I,-4)); % in fact, when the V.Bitdepth =16 in a hamamatsu multitif,the bitdepth=12
                else
                    I = uint8(I);
                end
            end
            if complement
                I = imcomplement(I);
            end
            I = I(pixels_range(1,1):pixels_range(1,2),...
                pixels_range(2,1):pixels_range(2,2));
            [I, ~] = image_normalization(I, background, minV, maxV, scaling);
            if ~isnan(boundary_threshold)
                % identifies only particles without holes
                [Boundaries{ff}, Centroids{ff}, Geometries{ff}]=...
                    bwboundaries_noholes(I, boundary_threshold,...
                    particle_threshold, minimum_area, false);
            end
            parfor_progress;
        end
        
    elseif strcmpi(ext, '.mat')
        fns = dir([path_name fn(1:end-5),'*.png']);
        parfor_progress(N);
        parfor ff = frames
            I = imread([path_name fns(ff).name]);
            I = I(pixels_range(1,1):pixels_range(1,2),...
                pixels_range(2,1):pixels_range(2,2));
            [I, ~] = image_normalization(I,background,minV,maxV, scaling);
            if ~isnan(boundary_threshold)
                % identifies only particles without holes
                [Boundaries{ff}, Centroids{ff}, Geometries{ff}]=...
                    bwboundaries_noholes(I, boundary_threshold,...
                    particle_threshold, minimum_area, false);
            end
            parfor_progress;
        end            
    end
end

Geometries = cellfun(@single,Geometries,'UniformOutput', 0);
Centroids = cellfun(@single,Centroids,'UniformOutput', 0);

fn = fullfile(path_name,[fn ext]);
save(outfile,'Boundaries', 'Centroids', 'Geometries', 'fn','-append')