function tracks = track_select(varargin)
% Finds tracks in cell array t (created by tracker.m) that are longer than
% threshold_length and creates a MX18XN matrix with all the good tracks. M
% is the maximum track length in all tracks, and N is the number of good
% tracks detected. Columns are:
%       1: track ID
%       2: frame number,
%       3: particle ID (within each frame)
%       4: time
%       5: x position of the centroid in microns
%       6: y position of the centroid in microns
%       7: area in microns^2
%       8: major axis in microns of the ellipse with the same normalized second
%           central moments as the region
%       9: minor axis in microns of the ellipse with the same normalized second
%           central moments as the region 
%       10: eccentricity
%       11: equivalent diameter
%       12: orientation in degrees
%       13: perimeter in microns
%       14: solidity
%       15: arc length in micronsm
%       16: minimum distance between poles of the particle in microns
%       17: is the width in microns 
%       18: is the radius of internal curvature in microns
%       19: is the R^2 of the fit to estimate cell curvature
%       20: is the d.f. of the fit to estimate cell curvature
%       21: is the RMSE of the fit to estimate cell curvature
%       22: is the amplitude of the sinusoidal fit (in microns)
%       23: is the frequency of the sinusoidal fit (in microns^-1)
%
% Saves data in a mat file and in a csv file (or and excel file, if you
% comment out indicated lines, with the same filename than the original
% video.
%
% input variables (optional):
%
%       fn = filename (optional)
%
%       threshold_length (optional) is the minimum length of the track to be
%           considered as such (default 2 seconds)
%
%       t (optional) is a cell array, in which every cell is a
%           MX2 matrix where M is the number of Frames that encodes a
%           particular track. The first column is the frame number, and the
%           second the particle ID in that given frame. N is the number of
%           frames of the track.
%
%       Centroids (optional) is a cell array of size (NumberofFramesX1), in
%           which each cell is a NX2 matrix, N being the number of
%           particles detected within the frame, and the 2 columns
%           represent position x and y in pixels of the centroids of the
%           particles.
% 
%       Geometries (optional) is a cell array of size (NumberofFramesX1),
%           in which each cell is a NX12 matrix, N being the number of
%           particles. Columns are area, MajorAxisLength, MinorAxisLength,
%           eccentricity, equivDiameter, orientation, perimeter and
%           solidity of each particle as defined in function regionprops,
%           and arc_length, minimum distance between poles, width and
%           interior radius of curvature.
%
%       Scaling (optional) is pixel/�m
%
%       FrameRate (optional) is the video frame rate in frames/s
%
%
% Requires: image processing toolbox, parfor_progress
%
%
% Ex.: [tracks, good_tracks, track_lengths] =...
%     track_select(fn, threshold_length, t, centroids, geometries, scaling, FrameRate)
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
%  phase-contrast microscope

%% input arguments
narginchk(0, 7)
if nargin==0 || isempty(varargin{1})
    [fn, fp] = uigetfile({'*.avi; *.mat'}, 'Choose an avi or mat file');
    fn = strcat(fp,fn);    
else
    fn = varargin{1};
end
[~,fn] = fileattrib(fn);
fn = fn.Name;

%% loads data
if nargin<2
    [pathstr,name] = fileparts(fn);
    matfile = fullfile(pathstr,[name,'.mat']);
    load(matfile, 't', 'Centroids', 'Geometries', 'scaling', '*Rate')
    threshold_length = 2;
elseif nargin<3
    [pathstr,name] = fileparts(fn);
    matfile = fullfile(pathstr,[name,'.mat']);
    load(matfile, 't', 'Centroids', 'Geometries', 'scaling', '*Rate')
    threshold_length = varargin{2};
else
    narginchk(7,7)
    threshold_length = varargin{2};
    t = varargin{3};
    Centroids = varargin{4};
    Geometries = varargin{5};
    scaling = varargin{6};
    FrameRate = varargin{7};
end

if exist('SamplingRate')
    FrameRate = SamplingRate;
end

%% removes short tracks
track_lengths = cellfun('size',t,1); % length of all detected tracks
good_tracks = track_lengths>=floor(FrameRate*threshold_length);
t = t(good_tracks);

tracks = single(nan(max(cellfun(@(x) size(x,1),t)), size(Geometries{find(~cellfun(@isempty,Geometries),1)},2)+6, size(t,2)));
for jj = 1:size(t,2)
    tracks(1:size(t{jj}),2:3,jj) = t{jj}; % frame number and particle ID
    for ii = 1:size(t{jj},1)
        tracks(ii,5:6,jj) = Centroids{t{jj}(ii,1)}(t{jj}(ii,2),:)/scaling; % position in microm
        tracks(ii,7,jj) = Geometries{t{jj}(ii,1)}(t{jj}(ii,2),1)/scaling^2; % area in microm^2
        tracks(ii,[8:9,11,13, 15:18, 21:22],jj) = Geometries{t{jj}(ii,1)}(t{jj}(ii,2), [2:3, 5, 7, 9:12, 15:16])/scaling; % geometry in microm
        tracks(ii,23,jj) = Geometries{t{jj}(ii,1)}(t{jj}(ii,2), 17)*scaling; % geometry in microm^-1
        tracks(ii,[10,12,14, 19, 20],jj) = Geometries{t{jj}(ii,1)}(t{jj}(ii,2), [4,6,8, 13, 14]); %  eccentricity, orientation and solidity
    end
end
clear ii jj t temp Geometries Centroids good_tracks track_lengths

%% splits tracks that are likely caused by crossing cells and turning cells
minimum_area = pi()*2*1;
divided_tracks = nan(size(tracks,1),size(tracks,2),0);
for jj = 1:size(tracks,3)
    kk = find(abs(diff(tracks(:,7,jj)))>minimum_area);% breaks the track when there is a change in area larger than 25%
%       kk = find(abs(diff(tracks(:,7,jj)))./tracks(1:end-1,7,jj)>0.25); % breaks the track when there is a change in area larger than 25%
    if ~isempty(kk)
        if kk(1)~=1
            kk = [1;kk];
        end
        k = find(~isnan(tracks(:,15,jj)),1,'last');
        if kk(end)~= k
           kk = [kk; k];
        end
        ff = flip(kk(2:end-1));
        for ii = 1:length(ff)
            if sum(~isnan(tracks(ff(ii):kk(end),5,jj)))>=floor(FrameRate*threshold_length)
              divided_tracks(:,:,end+1) = nan;
              divided_tracks(1:kk(end)-(ff(ii))+1,:,end) = tracks(ff(ii):kk(end),:,jj);
            end
            tracks(ff(ii):kk(end),:,jj) = nan;
        end
    end
end

tracks = cat(3,tracks, divided_tracks);

%% removes short tracks
good_tracks = squeeze(sum(~isnan(tracks(:,5,:))))>=floor(FrameRate*threshold_length);
tracks = tracks(:,:,good_tracks);

%% Save and export
[fp, fn] = fileparts(fn);
save(fullfile(fp,[fn '.mat']),...
    'tracks','FrameRate','threshold_length','scaling',...
    '-append')

%% xls file
% % comment this out for excel files
% n = size(tracks,3);
% for ii=1:n
%  xlswrite([fp, fn '.xls'],tracks(:,:,ii),sprintf('sheet%d',ii))
% end

%% csv file
% % comment this out for csv files
% fn = fullfile(fp,[fn '.csv']);
% fid = fopen(fn,'wt');
% fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
%     'track_ID','frame_number','particle_ID', 'time', 'xCentroid', 'yCentroid',...
%     'area', 'MajorAxis', 'MinorAxis', 'eccentricity', 'equivalent_diameter',...
%     'orientation', 'perimeter', 'solidity');
% fclose(fid);
% T = permute(tracks, [1,3,2]);
% T = reshape(T, size(tracks,1)*size(tracks,3), size(tracks,2));
% dlmwrite (fn, T, '-append');
