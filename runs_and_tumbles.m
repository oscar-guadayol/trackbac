function [run_lengths, run_speeds, run_complete,...
    run_angle_deviation, mean_run_orientation, run_path_frequency,...
    tumble_times, tumble_angles, tumble_angular_velocities,...
    speeds, angular_velocities, runs, tumbles, time, displacements,...
    angle_threshold, speed_threshold] =...
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
%           a change in direction to be considered a tumble. If empty, a
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
%           mean, median and max. Default is mean.
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

% [angle, displacement] = brownian_thresholds(Dr,Dt, win_length/FrameRate, 0.01);
% if isempty(angle_threshold) || isnan(angle_threshold)
%     angle_threshold = angle/(win_length/FrameRate); % apparent angular velocity threshold in radians per second
% end % if isempty(angle_threshold)
% if isempty(speed_threshold) || isnan(speed_threshold)
%     speed_threshold = displacement*1e6/(win_length/FrameRate); % apparent speed threshold in µmeters per second during time step
% end % if isempty(speed_threshold)

run_lengths = [];
run_speeds = [];
run_complete = [];
run_angle_deviation = [];
tumble_times = [];
tumble_angles = [];
tumble_angular_velocities = [];
mean_run_orientation = [];
minimum_run_length = minimum_run_length*FrameRate; % minimum run length in frames
displacements = [];
run_path_frequency = [];

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
    ftrack = smooth_track(track(:,5:6), win_length, method);
else
    ftrack = track(:,5:6);
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
    % regionprops uses the image with y axis reversed to calculate the
    % particle orientations. On top, the values range from -90 to 90.
    track(:,12) = mod((360-track(:,12))/180*pi,pi);
    
    % decomposes the angle
    y = sin(track(:,12));
    x = cos(track(:,12));
    
    reverse = all(sign([x y])~=sign([gradient(ftrack(:,1)),gradient(ftrack(:,2))]),2);
    x(reverse) = -x(reverse);
    y(reverse) = -y(reverse);
    
    reverse = find(any(sign([x y])~=sign([gradient(ftrack(:,1)),gradient(ftrack(:,2))]),2));
    for ii = 1:length(reverse)
        if reverse(ii) ==1
            jj = find(gradient(reverse)>1,1)+1;
        else
            jj = reverse(ii)-1;
        end
        
        if min(2*pi-abs(atan2(y(jj),x(jj))-atan2(y(reverse(ii)),x(reverse(ii)))),abs(atan2(y(jj),x(jj))-atan2(y(reverse(ii)),x(reverse(ii)))))>...
                min(2*pi-abs(atan2(y(jj),x(jj))-atan2(-y(reverse(ii)),-x(reverse(ii)))),abs(atan2(y(jj),x(jj))-atan2(-y(reverse(ii)),-x(reverse(ii)))))
            x(reverse(ii)) = -x(reverse(ii));
            y(reverse(ii)) = -y(reverse(ii));
        end
        
    end
    
    non_smoothed_orientations = atan2(y,x);
    non_smoothed_orientations(non_smoothed_orientations<0) = 2*pi+non_smoothed_orientations(non_smoothed_orientations<0);
    non_smoothed_orientations = unwrap(non_smoothed_orientations);
    
    % smooths the orientations by making a running vectorial average of
    % window win_length
    if win_length > 1
        y = smooth_track(y, win_length, method);
        x = smooth_track(x, win_length, method);
    end % if win_length > 1
else
    non_smoothed_orientations = atan2(y,x);
    non_smoothed_orientations(non_smoothed_orientations<0) = 2*pi+non_smoothed_orientations(non_smoothed_orientations<0);
    non_smoothed_orientations = unwrap(non_smoothed_orientations);
    
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

if rr(1)
    starts = [1;starts];
end
if rr(end)
    ends = [ends;length(rr)];
end
runs = [starts ends];
runs(diff(runs,[],2)<minimum_run_length &... % remove runs shorter than minimum run length
    runs(:,1)>1 &...        % and that are not at the beggining
    runs(:,2)<length(rr),:) = [];  % or the end of the track
tumbles = [runs(1:end-1,2)+1 runs(2:end,1)-1];

%% removes tumbles at the beggining and end of the track
tumbles(any(tumbles==1,2),:) = [];
tumbles(any(tumbles==length(rr),2),:) = [];

%% removes tumbles that are not flanked by full runs
% tumbles(tumbles(:,1) < min(runs(:,2)),:) = [];
% tumbles(tumbles(:,1) > max(runs(:,2)),:) = [];

%% statistics
tumble_times = diff(tumbles,[],2)/FrameRate; % run lengths in seconds
tumble_angles = nan(size(tumbles,1),1); % change in angle after a tumble
tumble_angular_velocities = nan(size(tumbles,1),1); % angular velocities during tumbles

for ii = 1:length(tumble_times)
    if tumbles(ii,1)>ceil(win_length/2)+1 &&...
            tumbles(ii,2)+1+floor(win_length/2)<length(time)
        angle0 = mean(orientations(tumbles(ii,1)-1-floor(win_length/2):tumbles(ii,1)-1));
        anglef = mean(orientations(tumbles(ii,2)+1:tumbles(ii,2)+1+floor(win_length/2)));
        tumble_angles(ii) = anglef-angle0;
    end
    tumble_angular_velocities(ii) = mean(angular_velocities(tumbles(ii,1):tumbles(ii,2)));
end

%% RUN ANGLE DEVIATIONS
if ~isempty(runs)
    run_lengths = (diff(runs,[],2)+1)/FrameRate; % run lengths in seconds
    run_complete = ~(any(runs==1,2) | any(runs==length(rr),2));
    run_speeds = nan(size(runs,1),4); % run mean speeds in µm/seconds
    
    run_angle_deviation = double.empty(0,50);
    mean_run_orientation = nan(size(runs,1),1);
    orientations = unwrap(orientations);
    displacements = nan(size(track,1), size(runs,1));
    for ii = 1:size(runs,1)
        displacements(1:runs(ii,2)-runs(ii,1)+1,ii) = ...
            ((ftrack(runs(ii,1),1)-ftrack(runs(ii,1):runs(ii,2),1)).^2+(ftrack(runs(ii,1),2)-ftrack(runs(ii,1):runs(ii,2),2)).^2);
        r_speeds = sqrt(sum(diff(ftrack(runs(ii,1):runs(ii,2),:),1,1).^2,2))*FrameRate;
        run_speeds(ii,1) = mean(r_speeds);
        run_speeds(ii,2) = std(r_speeds);
        min_speeds = min(r_speeds);
        max_speeds = max(r_speeds);
        mean_run_orientation(ii) = mean(orientations(runs(ii,1):runs(ii,2)));
        if isempty(min_speeds)
            run_speeds(ii,3) = nan;
            run_speeds(ii,4) = nan;
        else
            run_speeds(ii,3) = min_speeds;
            run_speeds(ii,4) = max_speeds;
        end
        
        run_angle_deviation = [run_angle_deviation;...
            orientations(runs(ii,1):runs(ii,2))-lagmatrix(orientations(runs(ii,1):runs(ii,2)),1:50)];
        
        run_path_frequency(ii)= nan;
        if runs(ii,2)-runs(ii,1)>8
            [ppx,f] = pwelch(detrend(non_smoothed_orientations(runs(ii,1):runs(ii,2))),diff(runs(ii,1:2))+1,[],[],FrameRate);
            [~, i] = max(ppx(2:end));
            run_path_frequency(ii)= f(i+1);
        end
        
    end
end
