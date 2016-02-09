function [run_lengths, run_speeds, run_angle_deviation, ...
    tumble_times, tumble_angles, tumble_angular_velocities,...
    speeds, angular_velocities, runs, tumbles, time] =...
    runs_and_tumbles(track, FrameRate, angle_threshold, speed_threshold,...
    minimum_run_length, a, b, method, win_length)

% Detects and analyses runs and tumbles within the input track, using the
% ellipsodial model for brownian motion to estimate angle and speed
% threhsolds to categorize segments of the track as runs or tumbles.
%
% Input variables:
%
%       track is a 2-column track with the (x,y) position of the particle
%           along the track.
%
%       FrameRate is the sampling rate of the video in frames per second.
%
%       angle_threshold (optional) is the threshold (in radians) for
%       aGamar1$!

%           change in direction to be considered a tumble. If empty, a
%           theoretical angle threshold is calculated using
%           brownian_thresholds.
%
%       speed_threshold (optional) is the threshold (in µm s^-1 ) for a
%       	section to be considered a run. If empty, a theoretical speed
%           threshold is calculated using brownian_thresholds.
%
%       minimum_run_length is the minimum length (in seconds) for a segment
%           of a track without significant  change in direction to be
%           considered as a run.
%
%       a and b are the semiaxis of the ellipsoid in meters.
%
%       method is the statistic used for the smoothing algorithm. It can be
%           mean, median and max. Default is method.
%
%       win_length is the window length (in frames) of the smoothing
%           algorithm.
%
% output variables:
%
%       run_lengths is vector with the length (in seconds) of all detected
%           runs.
%
%       run_speeds is a NX4 column where N is the number of complete runs
%           detected, 1st column is the mean speed, 2nd id the standard
%           deviation, 3rd is the minimum speed and 4th is the maximum
%           speed.
%
%       run_angle_deviation is a vector with the angle
%           deviations between frames separated by win_length of all runs
%           detected.
%
%       tumble_times is vector with the length (in seconds) of all
%           detected tumbles.
%
%       tumble_angles is vector with the angles (in radians) between
%           successive runs.
%
%       tumble_angular_velocities is a vector with the angular velocities
%           during tumbles.
%
%       speeds is a vector with the linear speeds of the track.
%
%       angular velocities is a vector with the angular velocities along
%           the track.
%
%       runs is a MX2 matrix where M is the number of runs detected in the
%           current track, and the first and second columns are the indexes
%           at which each particular run starts and ends respectevely.
%
%       tumbles is a MX2 matrix where M is the number of tumbles detected
%           in the current track, and the first and second columns are the
%           indexes at which each particular tumble starts and ends
%           respectevely.
%
%       time is a vector with the time in seconds of the track.
%
%
% Requires: image processing toolbox
%
% Called by: Runs_gui, Runs
%
% Ex.:
% [run_lengths, run_speeds, run_angle_deviation, ...
%     tumble_times, tumble_angles, tumble_angular_velocities,...
%     speeds, angular_velocities, runs, tumbles, time] =...
%     runs_and_tumbles(track, FrameRate, [], [],...
%     minimum_run_length, a, b, 'mean', win_length);
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

[angle, displacement] = brownian_thresholds(a/2,b/2, win_length/FrameRate, 0.01);
if isempty(angle_threshold) || isnan(angle_threshold)
    angle_threshold = angle/(win_length/FrameRate); % apparent angular velocity threshold in radians per second
end % if isempty(angle_threshold)
if isempty(speed_threshold) || isnan(speed_threshold)
    speed_threshold = displacement*1e6/win_length*FrameRate; % apparent speed threshold in µmeters per second during time step
end % if isempty(speed_threshold)

run_lengths = [];
run_speeds = [];
run_angle_deviation = [];
tumble_times = [];
tumble_angles = [];
tumble_angular_velocities = [];
minimum_run_length = minimum_run_length*FrameRate; % minimum run length in frames

%%
track = track(~isnan(track(:,5)),:);
time = (0:1/FrameRate:(size(track,1)-1)/FrameRate)';

if isempty(track)
    warning('No tracks')
    return
end

if isempty(method)
    method = 'mean';
end

%% Smoothing
if win_length > 1
    ftrack = smooth_track(track(:,5:6,:), win_length, method);
else
    ftrack = track(:,5:6,:);
end % if win_length > 1
x = gradient(ftrack(:,1),time);
y = gradient(ftrack(:,2),time);

%% Speeds
speeds = sqrt(x.^2+y.^2);

%% Angles
% the parameter used to estimate swimming angle depends on the aspect
% ratio. If aspect ratio is large, then it is assumed that the estimate of
% cell orientation from the particle detection script is robust, and it is
% used as angle estimator. Otherwise, the direction derived from the
% change in position of the centroid is used.

if a/b>3
        % this changes the -90 to 90 degrees output from the    % regionprops órientation'into a 0 to 180 degrees
    track(track(:,12)>0,12) = 180-track(track(:,12)>0,12);
    track(track(:,12)<0,12) = -track(track(:,12)<0,12);
    
    % this expands to 0 to 360 taking into account the direction of the
    % movement in the y axis. When the cell is running paralell to the x
    % axis, y oscillates between positive and negative, which produces
    % spurious reversals. To avoid this I use the y<-1 speed, rather than
    % the full precision value.
    track(y<-1,12)  = track(y<-1,12)+180;
    
    % decomposes the angle
    y = sin(track(:,12)/180*pi);
    x = cos(track(:,12)/180*pi);
    
    % smoothes the orientations by making a running vectorial average of
    % window win_length
    if win_length > 1
        y = smooth_track(y, win_length, method);
        x = smooth_track(x, win_length, method);
    end % if win_length > 1
    
end % if a/b>3

orientations = atan2(y,x);
orientations(orientations<0) = 2*pi+orientations(orientations<0);
angular_velocities = abs(gradient(unwrap(orientations),time));

%% runs detection
runs = [];
tumbles = [];
rr = angular_velocities<angle_threshold & speeds>speed_threshold;

starts = find((diff(rr))==1)+1;
ends = find((diff(rr))==-1);

if isempty(starts)||isempty(ends)
    return
end
if rr(1)
    starts = [1;starts];
end
if rr(end)
    ends = [ends;length(rr)];
end
runs = [starts ends];
runs(diff(runs,[],2)<minimum_run_length,:) = [];
tumbles = [runs(1:end-1,2)+1 runs(2:end,1)-1];

%% removes runs and tumble at the beggining and end of the track
runs(any(runs==1,2),:) = [];
runs(any(runs==length(rr),2),:) = [];
tumbles(any(tumbles==1,2),:) = [];
tumbles(any(tumbles==length(rr),2),:) = [];

%% statistics
tumble_times = diff(tumbles,[],2)/FrameRate; % run lengths in seconds
tumble_angles = nan(size(tumbles,1),1); % change in angle after a tumble
tumble_angular_velocities = nan(size(tumbles,1),1); % angular velocities during tumbles

for ii = 1:length(tumble_times)
    if tumbles(ii,1)>ceil(win_length/2) &&...
            tumbles(ii,2)+1+floor(win_length/2)<length(time)
        angle0 = mean(orientations(tumbles(ii,1)-1-floor(win_length/2):tumbles(ii,1)-1));
        anglef = mean(orientations(tumbles(ii,2)+1:tumbles(ii,2)+1+floor(win_length/2)));
        tumble_angles(ii) = anglef-angle0;
    end
    
    tumble_angular_velocities(ii) = ...
        atan2(...
        mean(sin(angular_velocities(tumbles(ii,1):tumbles(ii,2))/FrameRate)),...
        mean(cos(angular_velocities(tumbles(ii,1):tumbles(ii,2))/FrameRate))...
        )*FrameRate;
end

%% RUN ANGLE DEVIATIONS
if ~isempty(runs)
    run_lengths = diff(runs,[],2)/FrameRate; % run lengths in seconds
    run_speeds = nan(size(runs,1),4); % run mean speeds in µm/seconds
    run_angle_deviation = [];
    orientations = unwrap(orientations);    
    for ii = 1:size(runs,1)
        r_speeds = sqrt(sum(diff(ftrack(runs(ii,1):runs(ii,2),:)).^2,2))*FrameRate;
        run_speeds(ii,1) = mean(r_speeds);
        run_speeds(ii,2) = std(r_speeds);
        min_speeds = min(r_speeds);
        max_speeds = max(r_speeds);
        if isempty(min_speeds)
            run_speeds(ii,3) = nan;
            run_speeds(ii,4) = nan;
        else
            run_speeds(ii,3) = min_speeds;
            run_speeds(ii,4) = max_speeds;
        end
        time_step = win_length;
        track_deviation = orientations(runs(ii,1)+time_step:runs(ii,2)) -...
            orientations(runs(ii,1):runs(ii,2)-time_step);
        run_angle_deviation = [run_angle_deviation; track_deviation];
    end
end

%% Brownian thresholds
    function [angle, displacement] = brownian_thresholds(a,b, time_step, p_value)
        
        % Given an ellipsoid of revolution of semiaxis a and b, calculates
        % the maximum angle of rotational diffusivity around minor semiaxis
        % b and the maximum linear displacement (in m) along the major
        % semiaxis a expectable from brownian motion on a given time step.
        %
        % input variables:
        %       a and b are the semiaxis of the ellipsoid in meters
        %       time_step is the time between successive observations
        %       p_value (optional) is the probability of occurrance of a
        %           rotational angle wider than the threshold_angle.
        %           Default= 1-99.7300204/100 (the probability of 3
        %           standard deviations)
        %
        % output variables:
        %
        %       threshold angle around the minor aces (in radians)
        %       threshold dislacement along the major semiaxis (in meters)
        %
        % Oscar guadayol March 2015
        
        if nargin == 3 || isempty(p_value)
            p_value = 1-99.7300204/100;
        end
        
        [~,~, Dt, Dr]  = theoretical_friction_coefficients(a,b,b,33);
        
        
        % %% With flagellar stabilization, in Mitchell 2002 fashion
        % [ft_flagellum,fr_flagellum, ~, ~]  =...
        % theoretical_friction_coefficients(7.8*1e-6/2,0.2*1e-6/2,0.2*1e-6/2, 33);
        %
        % T = 33; % temperature
        % k_b = 1.38e-23; % Boltzmann constant (J/K)
        % K = T+273.15; % temperature in Kelvins
        %
        % Dt =  k_b*K./(ft + ft_flagellum);
        % Dr =  k_b*K./(fr + fr_flagellum);
        
        %%
        rms_angle = (2*Dr(2)*time_step).^(1/2); % standard deviation of the normal distribution of angles after a time=time_step
        angle = abs(norminv(p_value/2,0,rms_angle));
        
        rms_displacement = (2*Dt(1)*time_step).^(1/2); % standard deviation of the normal distribution of displacements after a time=time_step
        displacement = norminv(1-p_value/2,0,rms_displacement);
        
    end
end
