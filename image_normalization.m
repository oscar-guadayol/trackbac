function I = image_normalization(I, background, minI, maxI, denoise)
% Removes background, normalizes image to a 0-255 and subtracts noise
%
% input variables:
%
%       I is a grayscale image frame.
%
%       background is the average frame that is substracted from each frame
%       to get rid of non moving particles and noise
%
%       maxV and minV are the maximum and minimum intensity values of the
%           video after substracting the background and adding the mean
%           limits of the de-backgrounded images for histogram equalization
% 
%       denoise is a logical variable. Id true, each frame is denoised.
%
% output variables:
%
%       I is the normalized grayscale image.
%
% Called by particle_detection.m, threshold_level.m
%
% Requires: image processing toolbox
%
%   ex.:  I = image_normalization(I, background, minI, maxI, true)
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
%  phase-contrast microscope


% Subtracts video background
I = double(I)-double(background) + min(double(background(:))); 

% Adjusts image intensity values
I = (I-double(minI))*255/double(maxI-minI);

% Removes noise
if denoise
    SE = strel('disk',4,4);
    noise = imerode(imdilate(I,SE),SE);
    I = I-noise;
    I = (I-min(I(:)))*255/(max(I(:))-min(I(:)));
end
I = uint8(I);

