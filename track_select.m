function tracks = track_select(varargin)
% Finds tracks in cell array t (created by tracker.m) that are longer than
% threshold_length and creates a MX16XN matrix with all the good tracks. M
% is the maximum track length in all tracks, and N is the number of good
% tracks detected. First column is the track ID, 2nd is the frame number,
% 3rd is the particle ID (within each frame), 4th is time, 5th and 6th are
% x and y position of the centroid respectively in µm, 7th is the area in
% µm^2, 8th and 9th are the major and minor axes in µm of the ellipse that
% has the same normalized second central moments as the region, 10th the
% eccentricity, 11th is the equivalent diameter, 12th the orientation in
% degrees, 13th the perimeter in µm, 14th the solidity, 15 is the arc
% length and 16th is the minimum distance between poles of the particle
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
%       t (optional) is a cell array, in which every cell is a NumberofFramesX2 matrix
%           that ncodes a particular track. The first column is the frame
%           number, and the second the particle ID in that given frame. N
%           is the number of frames of the track.
%
%       Centroids (optional) is a cell array of size (NumberofFramesX1), in which each
%          cell is a NX2 matrix, N being the number of particles detected
%          within the frame, and the 2 columns represent position x and y
%          in pixels of the centroids of the particles.
%
%       Geometries (optional) is a cell array of size (NumberofFramesX1), in which
%           each cell is a NX8 matrix, N being the number of particles.
%           Columns are area, MajorAxisLength, MinorAxisLength,
%           eccentricity, equivDiameter, orientation, perimeter and
%           solidity of each particle as defined in function regionprops.
%
%       Scaling (optional) is pixel/µm
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
    load(matfile, 't', 'Centroids', 'Geometries', 'scaling', 'FrameRate')
    threshold_length = 2;
elseif nargin<3
    [pathstr,name] = fileparts(fn);
    matfile = fullfile(pathstr,[name,'.mat']);
    load(matfile, 't', 'Centroids', 'Geometries', 'scaling', 'FrameRate')
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

tracks = nan(size(Centroids,1), size(Geometries{find(~cellfun(@isempty,Geometries),1)},2)+6, size(t,2));
for jj = 1:size(t,2)
    tracks(t{jj}(:,1),2:3,jj) = t{jj}; % frame number and particle ID
    for ii = 1:size(t{jj},1)        
        tracks(t{jj}(ii,1),5:6,jj) = Centroids{t{jj}(ii,1)}(t{jj}(ii,2),:); % position
        tracks(t{jj}(ii,1),7:end,jj) = Geometries{t{jj}(ii,1)}(t{jj}(ii,2),:); % position
    end
    tracks(~isnan(tracks(:,2,jj)),1,jj) = jj;
    tracks(:,4,jj) = 0:1/FrameRate:(size(tracks,1)-1)/FrameRate;
end
clear ii jj

%% converts data from pixels to µm
tracks(:,[5:6,8:9,11,13, 15:17],:) = tracks(:,[5:6,8:9,11,13, 15:17],:)/scaling;
tracks(:,7,:) = tracks(:,7,:)/scaling^2;

%% removes short tracks
good_tracks = squeeze(sum(~isnan(tracks(:,5,:))))>=FrameRate*threshold_length;
tracks = tracks(:,:,good_tracks);

%% splits tracks that are likely caused by crossing cells and turning cells
for jj = 1:size(tracks,3)
    kk = find(abs(diff(tracks(~isnan(tracks(:,7,jj)),7,jj)))>4);
    if ~isempty(kk)
        if kk(1)~=1
            kk = [1;kk];
        end
        k = find(~isnan(tracks(:,15,jj)),1,'last');
        if kk(end)~= k
           kk = [kk; k];
        end
        for ii = flip(2:length(kk))
            tracks(:,:,end+1) = nan;
            tracks(1:k-kk(ii)+1,:,end) = tracks(kk(ii):k,:,jj);
            tracks(kk(ii):k,:,jj) = nan;
        end
    end
end

%% removes short tracks
good_tracks = squeeze(sum(~isnan(tracks(:,5,:))))>=FrameRate*threshold_length;
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
