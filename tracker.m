function t = tracker(varargin)

% This function automatically tracks multiple bacteria using 'Boundaries'
%   Cell array output from particles_threshold
%
% input variables (optional):
%
%       fn = filename (optional)
%
%       frames is a vector with frames to be used by the script. (optional)
%
%       Boundaries is a cell array of size (NumberofFramesX1) in which each
%           cell is another array of 2-column matrices, one for each
%           particle in the frame, with x,y points of the boundary of the
%           particle. (optional)
%
% output variables:
%
%       t is a cell array, in which every cell is a NumberofFramesX2 matrix
%       that encodes a particular track. The first column is the frame
%       number, and the second the particle ID in that given frame. N is
%       the number of frames of the track.
%
%
% Requires: image processing toolbox
% Follows: particle_detection
%
% Calls: video_background, threshold_level, preprocess,
%   bwboundaries_noholes, Progress
%
% Ex.: t = tracker(fn);
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
narginchk(0, 3)
if nargin==0 || isempty(varargin{1})
    [fn, fp] = uigetfile({'*.avi; *.mat; *.tif'},...
        'Choose an avi, multitiff, or mat file');
    fn = strcat(fp,fn);
else
    fn = varargin{1};
end
[~,fn] = fileattrib(fn);
fn = fn.Name;

if nargin<2
    %% loads particle detection data
    [pathstr,name] = fileparts(fn);
    matfile = fullfile(pathstr,[name,'.mat']);
    if ~exist(matfile,'file')
        error('No particle data. Please run first particle detection script')
    end
    load(matfile, 'frames', 'Boundaries')
else
    narginchk(3,3)
    frames = varargin{2};
    Boundaries = varargin{3};
end

%% analysis
mm = 0;
dummy = Boundaries;
t{1} = [];
for ff = frames(1):frames(end-1) % frame
    for pp = 1:length(Boundaries{ff}) % particle
        if isempty(dummy{ff}{pp})
            continue
        end
        mm = mm+1;
        e = pp;
        t{mm} = [ff,pp]; % [frame, particle]
        for kk = ff:length(frames)-1 % frame
            %% using for loop
%             E = false(length(Boundaries{kk}),1);
%             for bb = 1:length(Boundaries{kk+1})
%                 E(bb) = any([inpolygon(Boundaries{kk+1}{bb}(:,1),Boundaries{kk+1}{bb}(:,2),Boundaries{kk}{e}(:,1),Boundaries{kk}{e}(:,2));...
%                     inpolygon(Boundaries{kk}{e}(:,1),Boundaries{kk}{e}(:,2),Boundaries{kk+1}{bb}(:,1),Boundaries{kk+1}{bb}(:,2))]);
%             end
%             e = find(E);            
            
            %% using cellfun  
            if isempty(Boundaries{kk+1})
                break
            end
            C = cellfun(@(x)...
                any(inpolygon(x(:,1),x(:,2),Boundaries{kk}{e}(:,1),Boundaries{kk}{e}(:,2))),...
                Boundaries{kk+1},'UniformOutput',false);
            D = cellfun(@(x)...
                any(inpolygon(Boundaries{kk}{e}(:,1),Boundaries{kk}{e}(:,2),x(:,1),x(:,2))),...
                Boundaries{kk+1},'UniformOutput',false);
            
            % alternatives
            %             C = cellfun(@(x)...
            %                 overlap(x,Boundaries{kk}{e},h,w),...
            %                 Boundaries{kk+1},'UniformOutput',false);
            %             D = cellfun(@(x)...
            %                 overlap(x,Boundaries{kk+1}{e},h,w),...
            %                 Boundaries{kk},'UniformOutput',false);
            %             C = cellfun(@(x) any(inpoly(x,Boundaries{kk}{e})),...
            %                 Boundaries{kk+1},'UniformOutput',false);
            %             D = cellfun(@(x) any(inpoly(Boundaries{kk}{e},x)),...
            %                 Boundaries{kk+1},'UniformOutput',false);
            
            e = find(any(cell2mat([C,D]),2));
            
            if numel(e)==1
                t{mm} = [t{mm};kk+1,e]; % first column of each track is the frame, second is the particle
                dummy{kk+1}{e} = []; % to remove particles that cross each
                % and to avoid retracking already tracked particles
            elseif numel(e)>1
                for ll = 1:numel(e)
                    dummy{kk+1}{e(ll)} = []; % to remove particles that cross each
                    % and to avoid retracking already tracked particles
                end
                break
            else
                break
            end
        end
    end
end


%% 3D tracker. Work in progress
% from http://www.mathworks.com/matlabcentral/answers/85180-multi-dimensional-version-of-bwboundaries
% I = read(V);
% BW = imsubtract(I,repmat(background,1,1,size(I,3)))>uint8(boundary_threshold.dark);
% BW=((BW-min(BW(:))).*255./(max(BW(:))-min(BW(:))));
% BW = uint8(BW);
% BW = BW>boundary_threshold.dark;
% CC=bwconncomp(BW);
% Images = regionprops(CC,'Image','Centroid');
% t = Images.Image;

%% saving
save(matfile,'t','-append')