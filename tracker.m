function t = tracker(varargin)

% This function automatically tracks multiple bacteria using 'Boundaries'
%   Cell array output from particles_threshold
%
% input variables (optional):
%
%       fn = filename (optional)
%
%       frames, Boundaries, Centroids, and scaling are all
%       outputs of particle_detection.m (optional)
%
% output variables:
%
%       t is a cell array, in which every cell is a NumberofFramesX2 matrix
%           that encodes a particular track. The first column is the frame
%           number, and the second the particle ID in that given frame.
%
%
% Requires: image processing toolbox
% Follows: particle_detection
% Calls: Progress
%
%
% Requires: image processing toolbox, parfor_progress
%
% Calls: video_background, threshold_level, preprocess,
%   bwboundaries_noholes, Progress
%
% Ex.: t = tracker(fn);
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
narginchk(0, 5)
if nargin==0 || isempty(varargin{1})
    [fn, fp] = uigetfile({'*.avi; *.mat; *.tif'},...
        'Choose an avi, multitiff, or mat file');
    fn = strcat(fp,fn);
else
    fn = varargin{1};
end
[~,fn] = fileattrib(fn);
fn = fn.Name;

[pathstr,name] = fileparts(fn);
matfile = fullfile(pathstr,[name,'.mat']);
if nargin<2
    %% loads particle detection data
    if ~exist(matfile,'file')
        error('No particle data. Please run first particle detection script')
    end
    load(matfile, 'frames', 'Boundaries', 'Centroids', 'scaling')
else
    narginchk(5,5)
    frames = varargin{2};
    Boundaries = varargin{3};
    Centroids = varargin{4};
    scaling = varargin{5};
end

%% analysis
mm = 0;
dummy = Boundaries;
t{1} = [];
pr = Progress();
pr.push_bar('Track detection', 1, length(frames));
for ff = frames(1):frames(end-1) % frame
    pr.set_val(ff);
    for pp = 1:length(Boundaries{ff}) % particle
        
        if isempty(dummy{ff}{pp})
            continue
        end
        mm = mm+1;
        e = pp;
        t{mm} = [ff,pp]; % [frame, particle]
        for kk = ff:length(frames)-1 % frame
            %% using for loop
            % only looking at nearby particles in the next frame
            if isempty(Boundaries{kk+1})
                break
            end
            near = find(sqrt(sum((Centroids{kk+1}-repmat(Centroids{kk}(e,:),size(Centroids{kk+1},1),1)).^2,2))<=scaling*5); % finds particles in the following frame that are near the current particle (i,e, less than 5 microns away)
            E = false(length(near),1);
            for bb = 1:length(near)
                E(bb) = any([inpolygon(Boundaries{kk+1}{near(bb)}(:,1),Boundaries{kk+1}{near(bb)}(:,2),Boundaries{kk}{e}(:,1),Boundaries{kk}{e}(:,2));...
                    inpolygon(Boundaries{kk}{e}(:,1),Boundaries{kk}{e}(:,2),Boundaries{kk+1}{near(bb)}(:,1),Boundaries{kk+1}{near(bb)}(:,2))]);
            end             
            e = near(E);

            %% using cellfun  
            
%              % only looking at nearby particles in the next frame
%             near = find(sqrt(sum((Centroids{kk+1}-repmat(Centroids{kk}(e,:),length(Centroids{kk+1}),1)).^2,2))<=scaling*3);
%             C = cellfun(@(x)...
%                 any(inpolygon(x(:,1),x(:,2),Boundaries{kk}{e}(:,1),Boundaries{kk}{e}(:,2))),...
%                 Boundaries{kk+1}(near),'UniformOutput',false);
%             C = near(cell2mat(C));

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
%         toc
    end
end
pr.pop_bar();
t = cellfun(@uint16,t,'UniformOutput', 0);

%% saving
save(matfile,'t','-append')