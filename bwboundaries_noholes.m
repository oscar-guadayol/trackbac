function [boundaries, centroid, geometries] = bwboundaries_noholes(...
    I, boundary_threshold, particle_threshold, minimum_area, noholes)

% Detects particles in an 8 bit image that are above an intensity
% threshold.
%
% input variables:
%
%       I is a grayscale image matrix
%
%       boundary_threshold (optional) is the minimum pixel intensity in a
%           grayscale (0-255) of the boundaries of a particle.
%
%       particle_thresholds (optional) is the minimum pixel intensity in
%           a grayscale (0-255) that a particle needs to have to be
%           included in the analyses.
%
%       minimum_area is the minimum area in pixels, bacteria of 1µm radius
%           would have.
%
%       noholes is a logic variable. If true, bwboundaries searches for
%           object (parent and child) boundaries. This can provide better
%           performance. If false, bwboundaries searches for both object
%           and hole boundaries.
%
% output variables:
%
%       boundaries is a cell array of 2-column matrices, one for each
%           particle in the frame, with (x, y) positionsin pixels of
%           the boundary of the particle.
%
%       centroids is an NX2 matrix, N being the number of particles
%          detected within the frame, and the 2 columns represent (x,y)
%          position in pixels of the centroids of the particles.
%
%       geometries is an NX8 matrix, N being the number of particles.
%          Columns are area, MajorAxisLength, MinorAxisLength,
%          eccentricity, equivDiameter, orientation, perimeter and solidity
%          of each particle as defined in function regionprops.
%
% Called by particle_detection.m, threshold_level.m
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

centroid = [];
geometries = [];
boundaries = [];

%% applies a median filter with a minimum_area window to remove noise
win = round(sqrt(minimum_area)); % median filter window.
I = medfilt2(I,[win win]);

bw = im2bw(I,boundary_threshold/255);
bw = ~bw;

bw = imclearborder(bw); % removes particles on the edges
if ~any(bw(:))
    return
end

bw = bwareaopen(bw, minimum_area); % removes objects smaller than n pixels

if ~noholes
    [boundaries,L,n,A] = bwboundaries(bw);
    %% removes particles that have a hole inside
    donuts = nan(length(boundaries)-n,1);
    for ii = 1:length(donuts)
        donuts(ii) = find(A(ii+n,:));
    end
    
    % removes holes
    holes = (n+1:length(boundaries))';
    
    % marks particles in holes to be removed
    bad = [holes;donuts];
    
else
    [boundaries,L] = bwboundaries(bw,'noholes');
    bad = [];
end

for ii = 1:length(boundaries)
    if ~any(I(L==ii)<particle_threshold)
        bad = [bad; ii];
    end
end

boundaries(bad(:)) = [];
central_lines = bwboundaries(bwmorph(bw,'thin', Inf));
% central_lines(bad(:)) = [];
if isempty(boundaries)
    return
end

%% geometrical characterization
stats = regionprops(L, 'Centroid', 'Eccentricity', 'EquivDiameter', 'Area',...
    'MajorAxisLength', 'MinorAxisLength','Orientation', 'Perimeter', 'Solidity');
stats(bad) = [];

centroid = reshape([stats.Centroid],2,length(stats))';
geometries = nan(length(stats),11);
geometries(:,1) = [stats.Area]';
geometries(:,2) = [stats.MajorAxisLength]';
geometries(:,3) = [stats.MinorAxisLength]';
geometries(:,4) = [stats.Eccentricity]';
geometries(:,5) = [stats.EquivDiameter]';
geometries(:,6) = [stats.Orientation]';
geometries(:,7) = [stats.Perimeter]';
geometries(:,8) = [stats.Solidity]';

% removes geometrical parameters from phantom regions generated in
% bwboundaries_noholes
bad = [stats.Area]==0;
centroid(bad,:) = [];
geometries(bad,:) = [];

% [central_lines,~] = bwboundaries(bwmorph(bw,'thin', Inf),'noholes');
% central_lines(bad(:)) = [];
% extension = ceil(scaling);
extension = ceil(mean(geometries(:,2))/2);

for ii=1:length(boundaries)
    %% length
    % boundaries and arc_lengths are sometimes mismatched. this looks for the right arc_length for each boundary
    C = cellfun(@(x) inpolygon(x(1,1),x(1,2),boundaries{ii}(:,1),boundaries{ii}(:,2)),central_lines, 'UniformOutput',false);
    jj = find(cell2mat(C));
    
    % reorders arc length and removes repeated values
    [~,m,~] = unique(num2str(central_lines{jj}),'rows');
    M = central_lines{jj}(m,:);
    
    
    if size(M,1)<2
        geometries(ii,9) = range(boundaries{ii}(boundaries{ii}(:,2)==min(boundaries{ii}(:,2))+floor(range(boundaries{ii}(:,2))/2),1));
        geometries(ii,10) = range(boundaries{ii}(boundaries{ii}(:,2)==min(boundaries{ii}(:,2))+floor(range(boundaries{ii}(:,2))/2),1));
        geometries(ii,11) = range(boundaries{ii}(boundaries{ii}(:,1)==min(boundaries{ii}(:,1))+floor(range(boundaries{ii}(:,1))/2),2));
    else
        % expands the arc_length to the boundaries
        ll = shortest_path(M);
        len = find_arc_length(M(ll,:), boundaries{ii}, bw, 100);
        %         len = filtfilt((1/5)*ones(1,5), 1, len); % smoothes the central line
        geometries(ii,9) = sum(sqrt(sum(diff(len).^2,2))); % arc length of the particle
        geometries(ii,10) = sqrt(sum((len(1,:)-len(end,:)).^2)); % minimum distance between poles of the particle
        
        %% width
        if size(len,1)>2
            middle_point(1) = len(floor(size(len,1)/2),1)+len(floor(size(len,1)/2),1)-len(ceil(size(len,1)/2),1);
            middle_point(2) = len(floor(size(len,1)/2),2)+len(floor(size(len,1)/2),2)-len(ceil(size(len,1)/2),2);
            
            if mod((size(len,1)/2),2)==0
                middle_point(1) = len(size(len,1)/2,1)+(len(size(len,1)/2+1,1)-len(size(len,1)/2,1))/2;
                middle_point(2) = len(size(len,1)/2,2)+(len(size(len,1)/2+1,2)-len(size(len,1)/2,2))/2;
                alpha = atan2(...
                    diff(len(size(len,1)/2:size(len,1)/2+1,2)),...
                    diff(len(size(len,1)/2:size(len,1)/2+1,1)))...
                    +pi/2;
            else
                middle_point = len(ceil(size(len,1)/2),:);
                alpha = atan2(...
                    len(ceil(size(len,1)/2)-1,2) - len(ceil(size(len,1)/2)+1,2),...
                    len(ceil(size(len,1)/2)-1,1) - len(ceil(size(len,1)/2)+1,1))...
                    +pi/2;
            end
            h = -ceil(extension/2):1:ceil(extension/2);
            clear inter
            [inter(:,1), inter(:,2)] = intersections(h.*cos(alpha) + middle_point(1),...
                h.*sin(alpha) + middle_point(2),...
                boundaries{ii}(:,1),boundaries{ii}(:,2));
            inter = unique(round(inter*100)/100,'rows'); % removes duplicated intersections
            geometries(ii,11) = mean(pdist(inter));
        end
    end
end
end


%% Shortest path

function ll = shortest_path(M)
% finds the shortest path of a nX2 matrix if the contour is
% branching (for example, because two cells are crossing), it will
% choose one branch randomly

% finds the starting point as the point farthest away from anything
ll = nan(size(M,1),1);
D = squareform(pdist(M));
A = nan(length(D),2);
for ii = 1:length(D)
    A(ii,1) = sum(D(1:find(D(:,ii)==0),ii));
    A(ii,2) = sum(D(find(D(:,ii)==0):size(D,2),ii));
end
% ll(1) = find(sum(squareform(pdist(M)))==max(sum(squareform(pdist(M)))),1);
[ll(1),~] = find(A==max(A(:)),1);

for jj = 2:length(ll)
    N = M;
    N(ismember(1:length(ll),ll),:) = nan;
    MD = sum((N - repmat(M(ll(jj-1),:),length(N),1)).^2,2);
    md = find(MD == min(MD) & MD<=10);
    if isempty(md)
        ll(isnan(ll)) = [];
        return
    end
    if length(md)>1
        M(md(2:end),:) = nan;
    end
    ll(jj) = md(1);
end
end

function arc_length = find_arc_length(arc_length, boundary, bw, extension)

p = arc_length(end-1:end,end-1:end);
[row, column] = find_closest_edge(p, boundary, bw);
arc_length = [arc_length;row',column'];

p = flipud(arc_length(1:2,1:2));
[row, column] = find_closest_edge(p, boundary, bw);
arc_length = [row',column';arc_length];

    function [row,column] = find_closest_edge(p, boundary, bw)
        
        b = diff(p(:,2))/diff(p(:,1)); % slope
        a = p(1,2)-b*p(1,1); %intersect
        
        if isinf(b) % if the line is invariant in x (vertical line)
            if p(3)-p(4)<0
                y = min([p(4)+1,size(bw,2)]):min([p(4)+extension,size(bw,2)]);
                x = repmat(p(2),1,length(y));
                x = round(x);
                y = round(y);
                sub = sub2ind(size(bw),x,y);
                
                if bw(sub(1))
                    [row,column] = ind2sub(size(bw),sub(1:find(bw(sub)==0,1)-1));
                else
                    row = [];
                    column = [];
                end
                
                if ~isempty(row) &&...
                        ~any(ismember(sub2ind(size(bw), boundary(:,1), boundary(:,2)), sub(1:find(bw(sub)==0,1)-1)))
                    [r,c] = ind2sub(size(bw), sub(find(bw(sub)==0,1,'last'):find(bw(sub)==0,1,'last')+1));
                    row = [row, diff(r)/2+r(1)];
                    column = [column, diff(c)/2+c(1)];
                end
            elseif p(3)-p(4)>0
                y = max(p(4)-extension,1):max(p(4)-1,1);
                x = repmat(p(2),1,length(y));
                x = round(x);
                y = round(y);
                sub = sub2ind(size(bw),x,y);
                
                if bw(sub(end))
                    [row,column] = ind2sub(size(bw),sub(find(bw(sub)==0,1,'last')+1:find(bw(sub),1,'last')));
                else
                    row = [];
                    column = [];
                end
            end
        else
            if p(1)-p(2)<0
                x = p(2)+1:p(2)+extension;
                y = a + x*b;
                x = round(x);
                y = round(y);
                
                %% cuts x or y to avoid out of frame points
                if min(y)<1
                    x = x(y>=1);
                    y = y(y>=1);
                end
                if min(x)<1
                    y = y(x>=1);
                    x = x(x>=1);
                end
                if max(y)>size(bw,2)
                    x = x(y<=size(bw,2));
                    y = y(y<=size(bw,2));
                end
                if max(x)>size(bw,1)
                    y = y(x<=size(bw,1));
                    x = x(x<=size(bw,1));
                end
                
                %%
                sub = sub2ind(size(bw), x, y);
                if bw(sub(1))
                    [row,column] = ind2sub(size(bw),sub(find(bw(sub),1,'first'):find(bw(sub)==0,1)-1));
                else
                    row = [];
                    column = [];
                end
                
                if ~isempty(row) &&...
                        ~any(ismember(sub2ind(size(bw), boundary(:,1), boundary(:,2)), sub(1:find(bw(sub)==0,1)-1)))
                    [r,c] = ind2sub(size(bw), sub(find(bw(sub)==0,1)-1:find(bw(sub)==0,1)));
                    row = [row, diff(r)/2+r(1)];
                    column = [column, diff(c)/2+c(1)];
                end
                if isempty(row) &&...
                        ~any(ismember(sub2ind(size(bw), boundary(:,1), boundary(:,2)), sub2ind(size(bw),p(2),p(4))))
                    [r,c] = ind2sub(size(bw), sub(find(bw(sub)==0,1)));
                    row = (p(2)-r)/2+r;
                    column = (p(4)-c)/2+c;
                end
                
            elseif p(1)-p(2)>0
                x = p(2)-extension:p(2)-1;
                y = a + x*b;
                x = round(x);
                y = round(y);
                
                %% cuts x or y to avoid out of frame points
                if min(y)<1
                    x = x(y>=1);
                    y = y(y>=1);
                end
                if min(x)<1
                    y = y(x>=1);
                    x = x(x>=1);
                end
                if max(y)>size(bw,2)
                    x = x(y<=size(bw,2));
                    y = y(y<=size(bw,2));
                end
                if max(x)>size(bw,1)
                    y = y(x<=size(bw,1));
                    x = x(x<=size(bw,1));
                end
                
                %%
                sub = sub2ind(size(bw),x,y);
                if bw(sub(end))
                    [row,column] = ind2sub(size(bw),sub(find(bw(sub)==0,1,'last')+1:find(bw(sub),1,'last')));
                else
                    row = [];
                    column = [];
                end
                
                if ~isempty(row) &&...
                        ~any(ismember(sub2ind(size(bw), boundary(:,1), boundary(:,2)), sub(find(bw(sub)==0,1,'last')+1:find(bw(sub),1,'last'))))
                    [r,c] = ind2sub(size(bw), sub(find(bw(sub)==0,1,'last'):find(bw(sub)==0,1,'last')+1));
                    row = [diff(r)/2+r(1),row];
                    column = [diff(c)/2+c(1), column];
                end
                if isempty(row) &&...
                        ~any(ismember(sub2ind(size(bw), boundary(:,1), boundary(:,2)), sub2ind(size(bw),p(2),p(4))))
                    [r,c] = ind2sub(size(bw), sub(find(bw(sub)==0,1,'last')));
                    row = (p(2)-r)/2+r;
                    column = (p(4)-c)/2+c;
                end
            end
        end
    end
end
