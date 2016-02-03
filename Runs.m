function [Run_lengths, Run_speeds, Run_angle_deviations,...
    Tumble_times, Tumble_angles, Tumble_angular_velocities,...
    Cell_area, Cell_length, Cell_width, Cell_eccentricity,...
    Cell_equi_diam, Cell_perimeter, Cell_solidity,...
    StartTimes, fig] = ...
    Runs(tracks, fns, angle_threshold, speed_threshold, win_length,...
    minimum_run_length, method, FrameRate, plotting, figure_export)

% This function detects and analyses the runs and tumbles in all the input
% tracks and compiles all the pertinent statistics.
%
% Input variables:
%
%       tracks is the matrix output from track_select or merge_tracks
%
%       fns is a full path filename or a cell array in which elements are
%           the mat files full paths.
%
%       angle_threshold (optional) is the threshold (in radians) for a
%           change in direction to be considered a tumble. If empty,
%           runs_and_tumbles will estimate a theoretical threshold based on
%           the ellipsodial model from Brownian motion.
%
%       speed_threshold (optional) is the threshold (in µm s^-1 ) for a
%       	section to be considered a run. If empty, runs_and_tumbles will
%           estimate a theoretical threshold based on the ellipsodial
%           model from Brownian motion.
%
%       win_length is a scalar with the length (in frames) of the moving
%           average window applied to the x y positions of the track. 
%
%       minimum_run_length is the minimum length (in seconds) for a segment
%           of a track without significant  change in direction to be
%           considered as a run.
%
%       method is the statistic used for the smoothing algorithm. It can be
%           'mean', 'median' and 'max'. Default is 'mean'.
%
%       FrameRate is a scalar with the frame rate of the video in frames
%           per second.
%
%       plotting is a logical value. If true it will generate a multiple
%           axis plot showing the statistical distributions of of runs and
%           tumbles. Default is false.
%
%       figure_export is a logical value. If true a dialog will open to
%           save the figure as eps. 
%
%
% Output variables:
%       Run_lengths is an MX2 matrix, where M is the number of complete runs.
%           The first column is the track ID, and the second column is the
%           length of the run in seconds.
%
%       Run_speeds is an MX5 matrix, where M is the number of complete runs.
%           The first column is the track ID, the 2nd column is the mean
%           speed, 3rd is the standard deviation, 4th is the minimum speed
%           and 5th is the maximum speed.
%
%       Run_angle_deviations is an MX2 matrix, where the first column is
%           the track ID and the 2nd column is the deviations between
%           frames separated by win_length of all runs detected.
%
%       Tumble_times is vector with the length (in seconds) of all
%           detected tumbles.
%
%       Tumble_angles is vector with the angles (in radians) between
%           successive runs.
%
%       Tumble_angular_velocities is a vector with the angular velocities
%           during tumbles.
%
%       Cell_area is an MX2 matrix, where M is the number of complete runs.
%           The first column is mean area of the cell, the 2nd column is
%           the standard deviation in µm^2.
%
%       Cell_length is an MX2 matrix, where M is the number of complete runs.
%           The first column is mean length of the cell, the 2nd column is
%           the standard deviation in µm.
%
%       Cell_width is an MX2 matrix, where M is the number of complete runs.
%           The first column is mean width of the cell, the 2nd column is
%           the standard deviation in µm.
%
%       Cell_eccentricity is an MX2 matrix, where M is the number of
%           complete runs. The first column is mean excentricity of the
%           cell, the 2nd column is the standard deviation. See regionprops
%           for details.
%
%       Cell_equi_diam is an MX2 matrix, where M is the number of complete
%           runs. The first column is mean equivalent diameter of the cell,
%           the 2nd column is the standard deviation.  See regionprops for
%           details.
%
%       Cell_perimeter is an MX2 matrix, where M is the number of complete
%           runs. The first column is mean perimeter of the cell, the 2nd
%           column is the standard deviation. See regionprops for details.
%
%       Cell_solidity is an MX2 matrix, where M is the number of complete
%           runs. The first column is mean solidity of the cell, the 2nd
%           column is the standard deviation. See regionprops for details.
%
%       StartTimes is a vector with the starting times of the video for
%           each particle.
%
%       fig is the handle for the output figure.
%
%
% Requires: image processing toolbox
%
% Calls: runs_and_tumbles, runs_and_tumble_figure
%
% Called by: Runs_gui
%
% Ex.:
%   [Run_lengths, Run_speeds, Run_angle_deviations, Run_cos_deviations,...
%     Tumble_times, Tumble_angles, Tumble_angular_velocities,...
%     Cell_area, Cell_length, Cell_width, Cell_eccentricity,...
%     Cell_equi_diam, Cell_perimeter, Cell_solidity,...
%     StartTimes, fig] = ...
%     Runs(tracks, fns, angle_threshold, speed_threshold, win_length,...
%     minimum_run_length, 'mean', FrameRate, true, false)
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
%  phase-contrast microscope.

Run_lengths = [];
Run_speeds = [];
Run_angle_deviations = [];
Tumble_times = [];
Tumble_angles = [];
Tumble_angular_velocities = [];
Cell_area =  [];
Cell_length =  [];
Cell_width = [];
Cell_eccentricity =  [];
Cell_equi_diam = [];
Cell_perimeter = [];
Cell_solidity =  [];
StartTimes =  [];
fig = [];

for tt = 1:size(tracks,3)
    dimensions_offset = 0.7670; % offset bewteen the dimensions obtained with phase contrast and those from electon microscope
    a = (nanmedian(tracks(:,15,tt))-dimensions_offset)*1e-6;
    % b = (nanmedian(tracks(:,17,tt))-dimensions_offset)*1e-6;
    b = 0.71*1e-6;

   [run_lengths, run_speeds, run_angle_deviation, ...
       tumble_times, tumble_angles, tumble_angular_velocities] = ...
        runs_and_tumbles(tracks(:,:,tt), FrameRate, angle_threshold,...
        speed_threshold, minimum_run_length, a, b,...
        method, win_length, fns{tracks(1,end-1,tt)}(1:end-3), false, false);
    
    if any([~isempty(run_lengths), ~isempty(tumble_times)])
        run_lengths = [run_lengths; nan(length(tumble_times)-length(run_lengths),1)];
        run_speeds = [run_speeds; nan(length(tumble_times)-size(run_speeds,1),4)];
        run_angle_deviation = [run_angle_deviation; nan(length(tumble_times)-size(run_angle_deviation,1),1)];
        tumble_times = [tumble_times; nan(length(run_lengths)-length(tumble_times),1)];
        tumble_angles= [tumble_angles; nan(length(run_lengths)-length(tumble_angles),1)];
        tumble_angular_velocities= [tumble_angular_velocities; nan(length(run_lengths)-length(tumble_angular_velocities),1)];
        Run_lengths = [Run_lengths; [ones(length(run_lengths), 1)* tt, run_lengths]];
        Run_speeds = [Run_speeds; [ones(length(run_lengths), 1)* tt, run_speeds]];
        Run_angle_deviations = [Run_angle_deviations; [ones(length(run_angle_deviation), 1)* tt, run_angle_deviation]];
        Tumble_times = [Tumble_times; tumble_times];
        Tumble_angles = [Tumble_angles; tumble_angles];
        Tumble_angular_velocities = [Tumble_angular_velocities; tumble_angular_velocities];        
        Cell_area = [Cell_area; repmat([nanmean(tracks(:,7,tt)), nanstd(tracks(:,7,tt))],length(run_lengths),1)];
        Cell_length = [Cell_length; repmat([nanmedian(tracks(:,15,tt))-dimensions_offset, nanstd(tracks(:,15,tt))],length(run_lengths),1)];
        Cell_width = [Cell_width; repmat([nanmedian(tracks(:,17,tt))-dimensions_offset, nanstd(tracks(:,17,tt))],length(run_lengths),1)];
        Cell_eccentricity = [Cell_eccentricity; repmat([nanmean(tracks(:,10,tt)), nanstd(tracks(:,10,tt))],length(run_lengths),1)];
        Cell_equi_diam = [Cell_equi_diam; repmat([nanmean(tracks(:,11,tt)), nanstd(tracks(:,11,tt))],length(run_lengths),1)];
        Cell_perimeter = [Cell_perimeter; repmat([nanmean(tracks(:,13,tt)), nanstd(tracks(:,13,tt))],length(run_lengths),1)];
        Cell_solidity = [Cell_solidity; repmat([nanmean(tracks(:,14,tt)), nanstd(tracks(:,14,tt))],length(run_lengths),1)];
        if size(tracks,2)==19
            StartTimes = [StartTimes; repmat(tracks(1,end,tt),length(run_lengths),1)];
        end
    end
end

%% figure
if ~isempty(Run_speeds) && plotting
    fig = runs_and_tumbles_figure(Run_speeds, Run_lengths, Tumble_times, Tumble_angles, FrameRate, win_length, figure_export);
end
