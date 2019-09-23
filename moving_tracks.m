function  [moving, percentage_moving] =...
    moving_tracks(tracks, FrameRate, win_length, merged, dimensions_offset, a, b)

% Discriminates between actively swimming cells and cells moving purely by
%   Brownian motion.
%
% input variables:
%
%       tracks is a matrix as output from track_select or merge_tracks.
%
%       FrameRate is the frame rate of the video.
%
%       win_length is a scalar with the length (in frames) of the moving
%           average window applied to the x y positions of the track.
%
% output variables:
%
%       moving is a vector of logical values with length equal to the
%           number of tracks. 1s are actively swimming particles, 0s are
%           Brownianly moving particles.
%
% Called by merge_tracks
%
% Calls: theoretical_friction_coefficients.
%
% Requires: image processing toolbox, intersections
%
% Ex.:[boundaries, centroid, geometries] = bwboundaries_noholes(...
%       I, boundary_threshold, particle_threshold, minimum_area, noholes)
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

%% theoretical diffusivities
if nargin<5
    dimensions_offset = 0;
end

if ~exist('a','var') || isempty('a')
    a = squeeze(nanmedian(tracks(:,15,:))-dimensions_offset)*1e-6; % median length
end
if ~exist('b','var') || isempty('b')
    b = ones(length(a),1)* 0.75*1e-6;
end
nframes = size(tracks,1);
ntracks = size(tracks,3);
times = squeeze(tracks(:,end,:));
frames = squeeze(tracks(:,2,:));

tracks = [tracks(:,5:6,:), sin(tracks(:,12,:)), cos(tracks(:,12,:))];
[~,~, Dt, ~]  = theoretical_friction_coefficients(a./2,b./2,b./2, 30);

%% displacements along major axis
moving = false(ntracks,1);
if win_length > 1
    x = binner([(1:nframes)', squeeze(tracks(:,1,:))],1,1,nframes,win_length,0,@nanmean,1);
    y = binner([(1:nframes)', squeeze(tracks(:,2,:))],1,1,nframes,win_length,0,@nanmean,1);
    ax = binner([(1:nframes)', squeeze(tracks(:,3,:))],1,1,nframes,win_length,0,@nanmean,1);
    ay = binner([(1:nframes)', squeeze(tracks(:,4,:))],1,1,nframes,win_length,0,@nanmean,1);
    clear tracks
    dx = diff(x);
    dy = diff(y);
    orientation = atan2(ax, ay);
    displacement = sqrt(dx.^2+dy.^2)*1e-6;
    displacement_along_major_axis = displacement.*cos(atan(dy./dx)-orientation(1:end-1,:));
    displacement_along_major_axis = displacement_along_major_axis.*sign(dx); % this accounts for backward movements, but it is a temporary fix. Need to do this properly with trigonometrics
    
    parfor tt = 1:ntracks
        %% Approach 1: displacements along the major axis are distributed normally with var = 2*D*t
        t = displacement_along_major_axis(~isnan(displacement_along_major_axis(:,tt)),tt);
        %         t = dtrack(~isnan(dtrack(:,1)),1)*1e-6;
        if numel(t)>=round(0.5*FrameRate/win_length)
            try
                moving(tt) = ttest(t, 0, 'Alpha', 0.001) ||...
                    vartest(t, 2*Dt(tt,1)*(win_length/FrameRate), 'Tail','right','alpha', 0.001);
            end
        end
        
        %% Approach 2: Mean square distances at time step equal 4*D*t
        % MSD = nanmean(displacements_along_majoraxis.^2); % Mean square distances at time step 1/FrameRate
        % moving = MSD./4.*FrameRate>10*Dt(:,1)';
        
        %% Approach 3: The first lag of the autocorrelation function is non significant if the particle is moving in brownian motion
        %     t = dx(~isnan(dx(1:end-1,1,tt)),1,tt);
        %  %     [~,p] = corrcoef(t(1:end-1), t(2:end),'alpha', 0.01);
        %     moving(tt) = lbqtest(t, 5);
        % moving = logical(moving);
    end
    
else
    dtracks = diff(tracks(:,5:6,:));
    orientations = tracks(:,12,:).*pi/180;
end

%% percentage moving in each frame
if ~merged % one single video file
    m = frames(:,moving);
    m = m(:);
    %     n = tracks(:,2,:);
    frames = frames(:);
    percentage_moving = ...
        mean(histcounts(m,1:nframes)./histcounts(frames,1:nframes));
else % file of merged videos
    StartTimes = unique(times(1,:)); % this identifies the different videos in each list of fns.
    
    % Number of tracks in each frame and video.
    t = times;
    t = t(:);
    for tt = 1:length(StartTimes)
        t(t==StartTimes(tt)) = tt;
    end
    n = frames;
    n = n(:);
    t(isnan(n)) = [];
    n(isnan(n)) = [];
    N = full(sparse(n,t,1));
    
    % Number of moving tracks in each frame and video.
    if any(moving)
        t = times(:,moving);
        t = t(:);
        for tt = 1:length(StartTimes);
            t(t==StartTimes(tt)) = tt;
        end
        
        m = frames(:,moving);
        m = m(:);
        ii = isnan(m) | isnan(t);
        t(ii) = [];
        m(ii) = [];
        M = full(sparse(m,t,1));
        M = [M; zeros(size(N,1)-size(M,1), size(M,2))];
        M = [M, zeros(size(M,1), size(N,2)-size(M,2))];
    else
        M = zeros(size(N));
    end
    
    % ratio of moving to nonmoving tracks in each frame and video
    R = M./N;
    
    % average ratio of moving to nonmoving tracks in each video
    percentage_moving = nanmean(R);
end
