function [pI, I] = image_normalization(I, background, minI, maxI, scaling)
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
%       scaling in pixels/µm
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

%% Subtracts video background
pI = double(I)-double(background) + min(double(background(:))); 

%% Adjusts image intensity values
pI = (pI-double(minI))*255/double(maxI-minI);

%% Bandpass filter to remove noise
highpass_window = round(scaling*2);
pI = bpass(imcomplement(uint8(pI)),1,highpass_window);
pI = (pI -min(pI(:)))*255/double(max(pI(:))-min(pI(:)));
pI = imcomplement(uint8(pI));

% adds a black frame around corresponding to what has been lost during
% noise removal.
pI(:,1:highpass_window) = 0;
pI(1:highpass_window,:) = 0;
pI(size(I,1)-highpass_window:end,:) = 0;
pI(:,size(I,2)-highpass_window:end) = 0;

I(:,1:highpass_window) = 0;
I(1:highpass_window,:) = 0;
I(size(I,1)-highpass_window:end,:) = 0;
I(:,size(I,2)-highpass_window:end) = 0;