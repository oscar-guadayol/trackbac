function ftrack = smooth_track(track, win_length, method)
% Smooths track with a running mean, median or maximum according to value
% of variable method.
%
% input variables:
%
%       track is
%       win_length is the window length (in frames) of the smoothing
%           algorithm.
%
%       method is the statistic used for the smoothing algorithm. It can be
%           mean, median and max. Default is method.
%
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
%          each cell is a NX8 matrix, N being the number of particles.
%          Columns are area, MajorAxisLength, MinorAxisLength,
%          eccentricity, equivDiameter, orientation, perimeter and solidity
%          of each particle as defined in function regionprops.
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
%  This file is part of trackbac, a collection of matlab scripts to
%  geometrically characterize and track swimming bacteria imaged by a
%  phase-contrast microscope.

method = char(method);
padded_track = padarray(track,win_length,'replicate');
ftrack = nan(size(padded_track));
if strcmp(method, 'mean') && win_length > 1
    ftrack = filtfilt(ones(1,win_length)/win_length,1,padded_track);
elseif strcmp(method, 'median') && win_length > 1
    ftrack = medfilt1(padded_track,win_length);
elseif strcmp(method, 'max') &&  win_length > 1
    for ii = win_length:length(padded_track)-win_length
        ftrack(ii,:) = nanmax(padded_track(ii-(win_length-1)/2:ii+(win_length-1)/2,:));
    end
end
ftrack = ftrack(win_length+1:end-win_length,:);
end