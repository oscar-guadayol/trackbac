function  [moving, percentage_moving] =...
    moving_tracks(tracks, FrameRate, win_length)

% Discriminates between actively swimming cells and cells moving purely by
%   Brownian motion.
%
% input variables:
%
%       tracks is a matrix as output from track_select aor merge_tracks.
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
dimensions_offset = 0.7670;
a = squeeze(nanmedian(tracks(:,15,:))-dimensions_offset)*1e-6; % median length
b = ones(length(a),1)* 0.71*1e-6;
Dt = nan(size(tracks,3), 3);
parfor tt=1:size(tracks,3)
    [~,~, Dt(tt,:), ~]  = theoretical_friction_coefficients(a(tt)/2,b(tt)/2,b(tt)/2, 33);
end

%% displacements along major axis
% win_length = 1;
ftracks = nan(ceil(size(tracks,1)/win_length), 2, size(tracks,3));
orientations = nan(ceil(size(tracks,1)/win_length), 1, size(tracks,3));


if win_length > 1
    for tt = 1:size(tracks,3)
        track =  binner(...
            [(1:size(tracks,1))', tracks(:,5:6,tt), sin(tracks(:,12,tt)) cos(tracks(:,12,tt))],...
            1,1,size(tracks,1),win_length,0,@nanmean,1);
        ftracks(:,:,tt) = track(:,1:2);
        orientations(:,:,tt) = atan2(track(:,3), track(:,4));
    end
else
    ftracks = tracks(:,5:6,:);
    orientations = tracks(:,12,:).*pi/180;
end

dtracks = diff(ftracks);
displacements = sqrt(dtracks(:,1,:).^2+dtracks(:,2,:).^2)*1e-6;
displacements_along_majoraxis = displacements.*cos(atan(dtracks(:,2,:)./dtracks(:,1,:))-orientations(1:end-1,1,:));
displacements_along_majoraxis = displacements_along_majoraxis.*sign(dtracks(:,1,:)); % this accounts for backward moevements, but it is a temporary fix. Need to do this porperlty with trigonometrics
displacements_along_majoraxis = squeeze(displacements_along_majoraxis);
moving = false(size(tracks,3),1);

%% Approach 1: displacements along the major axis are distributed normally with var = 2*D*t
for tt = 1:size(tracks,3)
    t = displacements_along_majoraxis(~isnan(displacements_along_majoraxis(:,tt)),tt);
    if numel(t)>=round(2*FrameRate/win_length)
        moving(tt) = ttest(t, 0, 'Alpha', 0.001) ||...
            vartest(t, 2*Dt(tt,1)*(win_length/FrameRate), 'Tail','right','alpha', 0.001);
    end
end

%% Approach 2: Mean square distances at time step equal 4*D*t
% MSD = nanmean(displacements_along_majoraxis.^2); % Mean square distances at time step 1/FrameRate
% moving = MSD./4.*FrameRate>10*Dt(:,1)';

%% Approach 3: The first lag of the autocorrelation function is non significant if the particle is moving in brownian motion
% for tt = 1:size(tracks,3)
%     t = dx(~isnan(dx(1:end-1,1,tt)),1,tt);
%  %     [~,p] = corrcoef(t(1:end-1), t(2:end),'alpha', 0.01);
%     moving(tt) = lbqtest(t, 5);
% end
% moving = logical(moving);

%% percentage moving in each frame
if size(tracks,2)==17; % on single video file
    m = tracks(:,2,moving);
    m = m(:);
    n = tracks(:,2,:);
    n = n(:);
    percentage_moving = ...
        mean(histcounts(m,1:size(tracks,1))./histcounts(n,1:size(tracks,1)));
    
elseif size(tracks,2)==19;
    StartTimes = unique(tracks(1,19,:)); % this identifies the different cideos in each list of fns.

    % Number of moving tracks in each frame and video.
    t = tracks(:,end,moving);
    t = t(:);
    for tt = 1:length(StartTimes);
        t(t==StartTimes(tt)) = tt;
    end
    
    m = tracks(:,2,moving);
    m = m(:);
    ii = isnan(m) | isnan(t);
    t(ii) = [];
    m(ii) = [];
    M = full(sparse(m,t,1));
    
    % Number of tracks in each frame and video.
    t = tracks(:,end,:);
    t = t(:);
    for tt = 1:length(StartTimes);
        t(t==StartTimes(tt)) = tt;
    end
    n = tracks(:,2,:);
    n = n(:);
    t(isnan(n))=[];
    n(isnan(n))=[];
    N = full(sparse(n,t,1));
    
    % ratio of moving to nonmoving tracks in each frame and video
    R = M./N;
    
    % average ratio of moving to nonmoving tracks in each video
    percentage_moving = nanmean(R);
end
