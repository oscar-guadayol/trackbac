function [tracks, moving, percentage_moving] = merge_tracks(fns)
% This function merges all the tracks from the list of mat files in fns,
% and adds a column to the 'tracks' variable with the filename number.
%
% input variables:
%
%       fns is a cell array in which e;ements are the mat files full paths.
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

dim = nan(3,length(fns)); %
for bb = 1:length(fns)
    B{bb} = load(fns{bb});
    if ~isfield(B{bb},'Boundaries')
        warning(['particle_detection script has not been run on ' fns{bb} '; skipping file'])
        continue
    end
    if ~isfield(B{bb},'tracks')
        warning(['track_select.m script has not been run on ' B{bb}.fn '; skipping file'])
        continue
    end
    dim(1,bb) = size(B{bb}.tracks,1);
    dim(3,bb) = size(B{bb}.tracks,3);
end

%%
% creates a merged tracks matrix, with the first dimension being the
% maximum track length encountered in all the fns, the 2nd dimension is the
% number of columns per track plus 1 for the filenumber, and 3rd is the
% total numebr of tracks in all the files together.
d=0;
tracks = nan(nanmax(dim(1,:)),19,nansum(dim(3,:)));
for tt = 1:length(fns)
    tracks(1:dim(1,tt), 1:17,d+1:d+dim(3,tt)) = B{tt}.tracks(:,1:17,:);
    tracks(1:dim(1,tt), 18  ,d+1:d+dim(3,tt)) =  tt; % *ones(dim(1,tt),1,dim(3,tt)); % adds a column with the filename number
    if isfield('B', 'StartTime')
        tracks(1:dim(1,tt), 19  ,d+1:d+dim(3,tt)) =  B{tt}.StartTime;  % adds a column with the time f file creation
    elseif isfield('B', 'ext') || ~isempty(regexp(B{tt}.ext, 'tif', 'once'))
        [fp, fn] = fileparts(char(fns(tt)));
        V = imfinfo([fp filesep fn '.tif']);
        StartTime = datenum(V(1).ImageDescription(33:52));
        tracks(1:dim(1,tt), 19  ,d+1:d+dim(3,tt)) =  StartTime;  % adds a column with the time f file creation
    end
    d = d + dim(3,tt);
end

load(fns{1},'FrameRate')

[moving, percentage_moving] = moving_tracks(tracks, FrameRate, 0);

%% saving
saving_dialog = questdlg('Save merged file?');

if strcmp(saving_dialog, 'Yes')
    % default folder and filename
    fn = double(char(fns));
    fn = char(fn(1,all(bsxfun(@eq, fn, fn(1,:)))));
    [pn,fn] = fileparts(fn);
    [fn,pn] = uiputfile('*.mat;*.csv','Save data', [pn, filesep, fn, '_merged_runs.mat']);
    fn = fullfile(pn, fn);
    save(fn, 'tracks','fns','moving', 'percentage_moving', '-v7.3')
end
