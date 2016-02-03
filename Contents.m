%  'tracback' is a toolbox for tracking and caracterizing bacterial tracks
%  that relies on the ellipsoidal model.
%  Copyright (C) 2016,  Oscar Guadayol <oscar_at_guadayol.cat>
%
% LICENSE:
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
% Parent functions:
%
%   particle_detection: detects and characterizes moving particles in a 
%       video.
%   tracker: tracks multiple bacteria using 'Boundaries' cell array
%       from particles_threshold.
%   track_select: select long tracks from cell
%       array t (created by tracker.m) and reorganizes them into a matrix.
%   merge_tracks: merges tracks from a list of mat files into a single
%       matrix.
%   Runs: detects and analyses the runs and tumbles in all the input tracks
%       and compiles all the pertinent statistics.
%   Runs_gui: allows to visualize individual tracks and their decomposition
%       in runs and tumbles, and process all the tracks within a list of
%       mat files, and saves the output.
%
% Child functions:
%
%   bwboundaries_noholes: detects and characterizes particles in an 8-bit
%       image.
%   threshold_level: interactive determination of intensity thresholds for
%       particle detection.
%   image_normalization: removes background and noise and renormalizes
%       image to a 0-255.
%   video_background: computes mean value of the whole video file.
%   moving_tracks: discriminates between actively swimming cells and cells
%       moving purely by Brownian motion.
%   theoretical_friction_coefficients: calculates translational and
%       rotational friction coefficients .
%   binner: bins a data array into uniform bin sizes using a user-defined
%       bin type 
%   runs_and_tumble_figure: creates a figure with the frequency distribution
%       of runs and tumbles in the input tracks.
%   smooth_track: smooths track with a running mean, median or maximum
%       according to value of variable method. 
%
% REFERENCES (sources):
% Guadayol et al 2016
