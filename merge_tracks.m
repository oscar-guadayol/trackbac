function [tracks, moving, percentage_moving] = merge_tracks(fns, win_length, dimensions_offset)
% This function merges all the tracks from the list of mat files in fns,
% and adds a column to the 'tracks' variable with the filename number.
%
% input variables:
%
%       fns is a cell array in which e;ements are the mat files full paths.
%
%       win_length is a scalar with the length (in frames) of the moving
%           average window applied to the x y positions of the track.
%           Default 0.1667 s
%
%       dimension_offset is the offset bewteen the dimensions obtained with
%           phase contrast microscopy and theoretical dimensions of
%           standard beads. Default is 0.7550 microns
%
% output variables:
%
%       tracks is a matrix as described in 'track_select' with an added
%           column with the index corresponding to the file names in fns.
%
% Called by Runs_gui
%
% Calls: moving_tracks
%
% Requires: image processing toolbox, intersections
%
% [tracks, moving, percentage_moving] = merge_tracks(fns, win_length, dimensions_offset)
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

%%
% loads the data form the different files into B structure, and retrieves
% the size of the track matrix in each file

if isempty(fns)
    [fns, pathname, ~] = uigetfile('multiselect','on',[pwd, filesep,  '*.mat']);
    fns = fullfile(pathname, fns);
end % if isempty(fns)

if ischar(fns)
    fns = {fns};
end

if nargin<2 || isempty(win_length)
    load(fns{1},'FrameRate')
    win_length = round(5/30*FrameRate);
end
if nargin<3 || isempty(dimensions_offset)
    dimensions_offset = 0.7550;
end
if length(fns)==1
    load(fns{1},'tracks','*Rate', 'ext');
    tracks(:, end+1, :) =  1;
    dim = size(tracks);
    if ~exist('StartTime','var')
        load(fns{1},'StartTime')
        tracks(:, end+1, :) =  StartTime;  % adds a column with the time of file creation
    elseif exist('ext','var') && ~isempty(regexp(ext, 'tif', 'once'))
        if strcmp(ext, '.tif')
            [fp, fn] = fileparts(char(fns(1)));
            V = imfinfo([fp '/' fn '.tif']);
            StartTime = datenum(V(1).ImageDescription(33:52));
            tracks(:, end, :) =  StartTime;  % adds a column with the time f file creation
        end
    else
        tracks(1:dim(1,ff), end  ,d+1:d+dim(3,ff)) =  ff;
    end
    if ~exist('FrameRate','var')
        FrameRate = SamplingRate;
    end % if ~exist('FrameRate', 'var')
    [moving, percentage_moving] = moving_tracks(tracks, FrameRate, win_length, 0, dimensions_offset);
    fn = cell2mat(fns);
    [pn,fn] = fileparts(fn);
    save([pn, filesep, fn, '_merged_tracks.mat'], 'tracks','fns','moving', 'percentage_moving', '-v7.3')
    return
end

dim = nan(3,length(fns)); %
for ff = 1:length(fns)
    m = matfile(char(fns(ff)));
    fileinfo=who(m);
    if ~any(ismember(fileinfo,'Boundaries'))
        warning(['particle_detection script has not been run on '...
            fns{ff} '; skipping file'])
        continue
    end
    if ~any(ismember(fileinfo,'tracks'))
        warning(['track_select.m script has not been run on '...
            fns{ff} '; skipping file'])
        continue
    end
    dim(:,ff) = size(m,'tracks')';
    
end
ext = m.ext;

%%
% creates a merged tracks matrix, with the first dimension being the
% maximum track length encountered in all the fns, the 2nd dimension is the
% number of columns per track plus 1 for the filenumber, and 3rd is the
% total numebr of tracks in all the files together.
d = 0;
tracks = nan(nanmax(dim(1,:)),dim(2,1)+2,nansum(dim(3,:)));
for ff = 1:length(fns)
    t = load(fns{ff},'tracks');
    tracks(1:dim(1,ff), 1:size(t.tracks,2),d+1:d+dim(3,ff)) = t.tracks;
    tracks(1:dim(1,ff), end-1  ,d+1:d+dim(3,ff)) =  ff; % adds a column with the filename number
    if ~exist('StartTime','var')
        load(fns{ff},'StartTime')
        tracks(1:dim(1,ff), end  ,d+1:d+dim(3,ff)) =  StartTime;  % adds a column with the time f file creation
    elseif exist('ext','var') && ~isempty(regexp(ext, 'tif', 'once'))
        if strcmp(ext, '.tif')
            [fp, fn] = fileparts(char(fns(ff)));
            V = imfinfo([fp filesep fn '.tif']);
            StartTime = datenum(V(1).ImageDescription(33:52));
            tracks(1:dim(1,ff), end  ,d+1:d+dim(3,ff)) =  StartTime;  % adds a column with the time f file creation
        end
    else
        tracks(1:dim(1,ff), end  ,d+1:d+dim(3,ff)) =  ff;
    end
    d = d + dim(3,ff);
end

load(fns{1},'*Rate')
if ~exist('FrameRate','var')
    FrameRate = SamplingRate;
end % if ~exist('FrameRate', 'var')
[moving, percentage_moving] = moving_tracks(tracks, FrameRate, win_length, 1, dimensions_offset);

%% saving
fn = double(char(fns));
fn = char(fn(1,all(bsxfun(@eq, fn, fn(1,:)))));
[pn,fn] = fileparts(fn);
save([pn, filesep, fn, '_merged_tracks.mat'], 'tracks','fns','moving', 'percentage_moving', '-v7.3')
