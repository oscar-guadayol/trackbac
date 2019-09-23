function [Run_lengths, Run_speeds, Run_complete,...
    Run_angle_deviations, Run_mean_orientations, Run_path_frequencies, ...
    Tumble_times, Tumble_angles, Tumble_angular_velocities,...
    Cell_area, Cell_length, Cell_width, Cell_eccentricity,...
    Cell_equi_diam, Cell_perimeter, Cell_solidity,...
    Cell_radius,Cell_radius_R2, Cell_radius_df,Cell_radius_RMSE,...
    Cell_amplitude, Cell_frequency,...
    StartTimes, fig] = ...
    Runs(tracks, cell_width, angle_threshold, speed_threshold, win_length,...
    minimum_run_length, method, FrameRate, dimensions_offset, plotting, figure_export)

% This function detects and analyses the runs and tumbles in all the input
% tracks and compiles all the pertinent statistics.
%
% Input variables:
%
%       tracks is the matrix output from track_select or merge_tracks
%
%       angle_threshold (optional) is the threshold (in radians) for a
%           change in direction to be considered a tumble. If empty,
%           runs_and_tumbles will estimate a theoretical threshold based on
%           the ellipsodial model from Brownian motion.
%
%       speed_threshold (optional) is the threshold (in Âµm s^-1 ) for a
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
% Requires: image processing toolbox, SW_viscosity from Thermophysical
%       properties of seawater toolbox (http://web.mit.edu/seawater/)
%
% Calls: runs_and_tumbles, runs_and_tumble_figure
%
% Called by: Runs_gui
%
% Ex.:
%   [Run_lengths, Run_speeds, Run_complete,...
%     Run_angle_deviations, Run_mean_orientations, Run_path_frequencies, ...
%     Tumble_times, Tumble_angles, Tumble_angular_velocities,...
%     Cell_area, Cell_length, Cell_width, Cell_eccentricity,...
%     Cell_equi_diam, Cell_perimeter, Cell_solidity,...
%     StartTimes, Displacements, fig] = ...
%     Runs(tracks, cell_width, angle_threshold, speed_threshold, win_length,...
%     minimum_run_length, method, FrameRate, dimensions_offset, true, false)
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

if isempty(cell_width)
    cell_width = 0.75*1e-6;
end
    
Run_lengths = cell(size(tracks,3),1);
Run_path_frequencies = cell(size(tracks,3),1);
Run_speeds = cell(size(tracks,3),1);
Run_complete = cell(size(tracks,3),1);
Run_angle_deviations = cell(size(tracks,3),1);
Run_mean_orientations = cell(size(tracks,3),1);
Tumble_times = cell(size(tracks,3),1);
Tumble_angles = cell(size(tracks,3),1);
Tumble_angular_velocities = cell(size(tracks,3),1);
Cell_area = cell(size(tracks,3),1);
Cell_length = cell(size(tracks,3),1);
Cell_width = cell(size(tracks,3),1);
Cell_eccentricity = cell(size(tracks,3),1);
Cell_equi_diam = cell(size(tracks,3),1);
Cell_perimeter = cell(size(tracks,3),1);
Cell_solidity = cell(size(tracks,3),1);
Cell_radius = cell(size(tracks,3),1);
Cell_radius_R2 = cell(size(tracks,3),1);
Cell_radius_df = cell(size(tracks,3),1);
Cell_radius_RMSE = cell(size(tracks,3),1);
Cell_amplitude = cell(size(tracks,3),1);
Cell_frequency = cell(size(tracks,3),1);

StartTimes = cell(size(tracks,3),1);
fig = [];

%% thresholds
A = (squeeze(nanmedian(tracks(:,15,:)))-dimensions_offset)*1e-6;
B = (squeeze(nanmedian(tracks(:,17,:)))-dimensions_offset)*1e-6;
[~,~, Dt, Dr]  = theoretical_friction_coefficients(A/2,repmat(cell_width/2,size(A,1),1),repmat(cell_width/2,size(A,1),1),30); 
p_value = 1-99.7300204/100;
time_step = win_length/FrameRate;
rms_angle = (2*Dr(:,2)*time_step).^(1/2); % standard deviation of the normal distribution of angles after a time=time_step
angle = abs(norminv(p_value/2,0,rms_angle));
rms_displacement = (2*Dt(:,1)*time_step).^(1/2); % standard deviation of the normal distribution of displacements after a time=time_step
displacement = norminv(1-p_value/2,0,rms_displacement);
if isempty(angle_threshold) || isnan(angle_threshold)
    angle_threshold = angle/(win_length/FrameRate); % apparent angular velocity threshold in radians per second
end % if isempty(angle_threshold)
if isempty(speed_threshold) || isnan(speed_threshold)
    speed_threshold = displacement*1e6/(win_length/FrameRate); % apparent speed threshold in µmeters per second during time step
end % if isempty(speed_threshold)

parfor tt = 1:size(tracks,3)
    track = tracks(:,:,tt);
    [Run_lengths{tt}, Run_speeds{tt}, Run_complete{tt},...
        Run_angle_deviations{tt}, Run_mean_orientations{tt}, Run_path_frequencies{tt},...
        Tumble_times{tt}, Tumble_angles{tt}, Tumble_angular_velocities{tt}, ...
        ~, ~, ~, ~, ~, ~,~,~] =...
        runs_and_tumbles(track, FrameRate, angle_threshold(tt),...
        speed_threshold(tt), minimum_run_length, A(tt), B(tt),...
        method, win_length);
    if any([~isempty(Run_lengths{tt}), ~isempty(Tumble_times{tt})])
        
        %% pads variables to equal lengths        
        Run_lengths{tt} = [Run_lengths{tt}; nan(length(Tumble_times{tt})-length(Run_lengths{tt}),1)];
        Run_path_frequencies{tt} = [Run_path_frequencies{tt}(:); nan(length(Tumble_times{tt})-length(Run_lengths{tt}),1)];
        Run_speeds{tt} = [Run_speeds{tt}; nan(length(Tumble_times{tt})-size(Run_speeds{tt},1),4)];
        Run_angle_deviations{tt} = [Run_angle_deviations{tt}; nan(length(Tumble_times{tt})-size(Run_angle_deviations{tt},1),1)];
        Run_angle_deviations{tt} = [ones(size(Run_angle_deviations{tt},1), 1)* tt, Run_angle_deviations{tt}];
        if isempty(Run_mean_orientations{tt})
            Run_mean_orientations{tt} = nan;
        end
        Tumble_times{tt} = [Tumble_times{tt}; nan(length(Run_lengths{tt})-length(Tumble_times{tt}),1)];
        Tumble_angles{tt} = [Tumble_angles{tt}; nan(length(Run_lengths{tt})-length(Tumble_angles{tt}),1)];
        Tumble_angular_velocities{tt} = [Tumble_angular_velocities{tt}; nan(length(Run_lengths{tt})-length(Tumble_angular_velocities{tt}),1)];
       
        %% add a column with the track number
        Run_lengths{tt} = [ones(length(Run_lengths{tt}), 1)* tt, Run_lengths{tt}];
        Run_complete{tt} = [ones(length(Run_complete{tt}), 1)* tt, Run_complete{tt}];
        Run_path_frequencies{tt} = [ones(length(Run_path_frequencies{tt}), 1)* tt, Run_path_frequencies{tt}];
        Run_speeds{tt} = [ones(size(Run_speeds{tt},1), 1)* tt, Run_speeds{tt}];
        Run_mean_orientations{tt} = [ones(length(Run_mean_orientations{tt}), 1)* tt, Run_mean_orientations{tt}];
        
        %% geometry of cells for each run/tumble
        Cell_area{tt} = repmat([nanmean(track(:,7)), nanstd(track(:,7))],size(Run_lengths{tt},1),1);
        Cell_length{tt} = repmat([nanmedian(track(:,15))-dimensions_offset, nanstd(track(:,15))],size(Run_lengths{tt},1),1);
        Cell_width{tt} = repmat([nanmedian(track(:,17))-dimensions_offset, nanstd(track(:,17))],size(Run_lengths{tt},1),1);
        Cell_eccentricity{tt} = repmat([nanmean(track(:,10)), nanstd(track(:,10))],size(Run_lengths{tt},1),1);
        Cell_equi_diam{tt} = repmat([nanmean(track(:,11)), nanstd(track(:,11))],size(Run_lengths{tt},1),1);
        Cell_perimeter{tt} = repmat([nanmean(track(:,13)), nanstd(track(:,13))],size(Run_lengths{tt},1),1);
        Cell_solidity{tt} = repmat([nanmean(track(:,14)), nanstd(track(:,14))],size(Run_lengths{tt},1),1);
        Cell_radius{tt} = repmat([nanmean(track(:,18)), nanstd(track(:,18))],size(Run_lengths{tt},1),1);
        Cell_radius_R2{tt} = repmat([nanmean(track(:,19)), nanstd(track(:,19))],size(Run_lengths{tt},1),1);
        Cell_radius_df{tt} = repmat([nanmean(track(:,20)), nanstd(track(:,20))],size(Run_lengths{tt},1),1);
        Cell_radius_RMSE{tt} = repmat([nanmean(track(:,21)), nanstd(track(:,21))],size(Run_lengths{tt},1),1);
        Cell_amplitude{tt} = repmat([nanmean(track(:,22)), nanstd(track(:,22))],size(Run_lengths{tt},1),1);
        Cell_frequency{tt} = repmat([nanmean(track(:,23)), nanstd(track(:,23))],size(Run_lengths{tt},1),1);

        StartTimes{tt} = repmat(track(1,end),size(Run_lengths{tt},1),1);
    end
end

Run_angle_deviations = vertcat(Run_angle_deviations{:});
Run_lengths = vertcat(Run_lengths{:});
Run_complete = vertcat(Run_complete{:});
Run_path_frequencies = vertcat(Run_path_frequencies{:});
Run_speeds = vertcat(Run_speeds{:});
Run_mean_orientations = vertcat(Run_mean_orientations{:});
Tumble_times = vertcat(Tumble_times{:});
Tumble_angles = vertcat(Tumble_angles{:});
Tumble_angular_velocities = vertcat(Tumble_angular_velocities{:});
Cell_area = vertcat(Cell_area{:});
Cell_length = vertcat(Cell_length{:});
Cell_width = vertcat(Cell_width{:});
Cell_eccentricity = vertcat(Cell_eccentricity{:});
Cell_equi_diam = vertcat(Cell_equi_diam{:});
Cell_perimeter = vertcat(Cell_perimeter{:});
Cell_solidity = vertcat(Cell_solidity{:});
Cell_radius = vertcat(Cell_radius{:});
Cell_radius_R2 = vertcat(Cell_radius_R2{:});
Cell_radius_df = vertcat(Cell_radius_df{:});
Cell_radius_RMSE = vertcat(Cell_radius_RMSE{:});
Cell_amplitude = vertcat(Cell_amplitude{:});
Cell_frequency = vertcat(Cell_frequency{:});
StartTimes = vertcat(StartTimes{:});

%% figure
if ~isempty(Run_speeds) && plotting
    fig = runs_and_tumbles_figure(Run_speeds, Run_lengths, Run_complete, Tumble_times, Tumble_angles, FrameRate, win_length, figure_export);
end
