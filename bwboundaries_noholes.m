function [boundaries, centroid, geometries] = bwboundaries_noholes(...
    I, boundary_threshold, particle_threshold, minimum_area, noholes, curvature, spirality)

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
%       minimum_area is the minimum area, in pixels, that a bacterium of
%           1 micron radius would have.
%
%       noholes is a logic variable. If true, bwboundaries searches for
%           object (parent and child) boundaries. This can provide better
%           performance. If false, bwboundaries searches for both object
%           and hole boundaries.
%
%       curvature is a logic variable. If true, bwboundaries calculates the
%           radius of curvature of each particle. Default is false.
%
%       spirality is a logic variable. If true, bwboundaries calculates the
%           amplitude and frequency of the cell by fitting a sinusoidal
%           function. Default is false.
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
%       geometries is an NX17 matrix, N being the number of particles.
%          Columns are area (1), MajorAxisLength (2), MinorAxisLength(3),
%          eccentricity(4), ebestquivDiameter(5), orientation(6),
%          perimeter(7) and solidity(8) of each particle as defined in
%          function regionprops, and arc_length(9), minimum distance
%          between poles(10), width(11), curvature radius(12:15) of the
%          best fit circle with goodness of fit parameters (R^2, degres of
%          freedom and RMSE), and amplitude(16) and frequency(17) of the
%          best fit sinusoid.
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
%  This file is part of trackbac v1.1, a collection of matlab scripts to
%  geometrically characterize and track swimming bacteria imaged by a
%  phase-contrast microscope

centroid = [];
geometries = [];
boundaries = [];

if ~exist('curvature') || isempty(curvature)
    curvature = false;
end
if ~exist('spirality') || isempty(spirality)
    spirality = false;
end

%% applies a median filter with a minimum_area window to remove noise

bw = imbinarize(I,boundary_threshold/255);
bw = ~bw;
bw = imclearborder(bw); % removes particles on the edges
if ~any(bw(:))
    return
end
bw = bwareaopen(bw, minimum_area); % removes objects smaller than n pixels

if ~noholes
    [boundaries, L] = bwboundaries(bw);
    bad = [];
    v = cell2mat(cellfun(@(v) v(1,:), boundaries, 'UniformOutput',false));  % matrix with the first values in each boundary
    for bb = 1:length(boundaries)
        A = find(inpolygon(v(:,1),v(:,2),boundaries{bb}(:,1),boundaries{bb}(:,2)));
        if numel(A)>1
            bad = [bad; A];
        end
    end
    
    % remove particles which minimum intensity is higher than the particle
    % threshold
    bad = [bad; setdiff((1:length(boundaries))', L(I<particle_threshold&L>0))];
    for ii = 1:length(boundaries)
        if ~any(I(L==ii)<particle_threshold)
            bad = [bad; ii];
        end
    end
    
    boundaries(bad(:)) = [];
    
    central_lines = bwboundaries(bwmorph(bw,'thin', Inf));
    v = cell2mat(cellfun(@(v) v(1,:), central_lines, 'UniformOutput',false)); % matrix with the first values in each central_lines
    C = cellfun(@(x) find(inpolygon(v(:,1),v(:,2),x(:,1),x(:,2))), boundaries, 'UniformOutput',false);
    central_lines = central_lines(cell2mat(C));
    
else
    [boundaries,L] = bwboundaries(bw,'noholes');
    bad = [];
    
    central_lines = bwboundaries(bwmorph(bw,'thin', Inf));
    v = cell2mat(cellfun(@(v) v(1,:), central_lines, 'UniformOutput',false)); % matrix with the first values in each central_lines
    C = cellfun(@(x) find(inpolygon(v(:,1),v(:,2),x(:,1),x(:,2)),1), boundaries, 'UniformOutput',false);
    central_lines = central_lines(cell2mat(C));
    
    % remove particles which minimum intensity is higher than the particle
    % threshold
    bad = [bad; setdiff((1:length(boundaries))', L(I<particle_threshold&L>0))];
    for ii = 1:length(boundaries)
        if ~any(I(L==ii)<particle_threshold)
            bad = [bad; ii];
        end
    end
    
    boundaries(bad(:)) = [];
    central_lines(bad(:)) = [];
end

if isempty(boundaries)
    return
end

%% geometrical characterization
stats = regionprops(L, 'Centroid', 'Eccentricity', 'EquivDiameter', 'Area',...
    'MajorAxisLength', 'MinorAxisLength','Orientation', 'Perimeter'); % remove solidity bevause it takes the most time in regionprops and I don't really use it. OG20170912
stats(bad) = [];

centroid = reshape([stats.Centroid],2,length(stats))';
geometries = nan(length(stats),17);
geometries(:,1) = [stats.Area]';
geometries(:,2) = [stats.MajorAxisLength]';
geometries(:,3) = [stats.MinorAxisLength]';
geometries(:,4) = [stats.Eccentricity]';
geometries(:,5) = [stats.EquivDiameter]';
geometries(:,6) = [stats.Orientation]';
geometries(:,7) = [stats.Perimeter]';

% removes geometrical parameters from phantom regions generated in
% bwboundaries_noholes
bad = [stats.Area]==0;
centroid(bad,:) = [];
geometries(bad,:) = [];

%%
extension = ceil(mean(geometries(:,2))/2);
h = -ceil(extension/2):1:ceil(extension/2);
for ii = 1:length(boundaries)
    
    % reorders arc length and removes repeated values
    [~,m,~] = unique(central_lines{ii},'rows');
    M = central_lines{ii}(m,:);
    
    if size(M,1)<2
        geometries(ii,9) = range(boundaries{ii}(boundaries{ii}(:,2)==min(boundaries{ii}(:,2))+floor(range(boundaries{ii}(:,2))/2),1));
        geometries(ii,10) = range(boundaries{ii}(boundaries{ii}(:,2)==min(boundaries{ii}(:,2))+floor(range(boundaries{ii}(:,2))/2),1));
        geometries(ii,11) = range(boundaries{ii}(boundaries{ii}(:,1)==min(boundaries{ii}(:,1))+floor(range(boundaries{ii}(:,1))/2),2));
        geometries(ii,12:15) = nan;
    else
        % expands the arc_length to the boundaries
        ll = shortest_path(M);
        len = find_arc_length(M(ll,:), boundaries{ii}, 100);
        len = [len(1,:); M(ll,:); len(end,:)];
        len([false; sum(diff(len).^2,2)==0],:) = []; % removes repeated consecutive pairs of values without resorting
        geometries(ii,9) = sum(sqrt(sum(diff(len).^2,2))); % arc length of the particle
        geometries(ii,10) = sqrt(sum((len(1,:)-len(end,:)).^2)); % minimum distance between poles of the particle
        
        %% width
        if size(len,1)>2
            clear inter ilen
            ilen = interp1([0; cumsum(sqrt(sum(diff(len).^2,2)))],len,0:0.1:(max(cumsum(sqrt(sum(diff(len).^2,2)))))','pchip');
            %             [~,middlepoint] = min(sum(abs(len-repmat(ilen(floor(size(ilen,1)/2),:),size(len,1),1)),2));
            [~,middlepoint] = min(sum((len(2:end-1,:)-repmat(ilen(floor(size(ilen,1)/2),:),size(len,1)-2,1)).^2,2));
            middlepoint = middlepoint+1;
            alpha = atan2(len(middlepoint-1,2) - len(middlepoint+1,2), len(middlepoint-1,1) - len(middlepoint+1,1))+pi/2;
            [inter(:,1), inter(:,2)] = intersections(h.*cos(alpha) + len(middlepoint,1),...
                h.*sin(alpha) + len(middlepoint,2),...
                boundaries{ii}(:,1),boundaries{ii}(:,2));
            
            geometries(ii,11) = mean(pdist(inter));
        end
        
        %% Curvature
        if curvature
            [radius, ~, ~, R2, df, RMSE] = radius_of_curvature(len(:,1), len(:,2));
            geometries(ii,12:15) = [radius-geometries(ii,11)/2,...
                R2, df, RMSE];
        else
            geometries(ii,12:15) = nan;
        end
        
        %% sinusoidal fit
        if spirality
            xy = central_lines{ii};
            xy = xy(1:ceil(size(xy,1)/2),:);
            s = warning('error', 'MATLAB:polyfit:RepeatedPointsOrRescale');
            try
                [amplitude, frequency, ~, ~] = sinusoidal(xy);
                geometries(ii,16) = amplitude;
                geometries(ii,17) = frequency;
            catch
                geometries(ii,16) = nan;
                geometries(ii,17) = nan;
            end
            warning(s)
        else
            geometries(ii,16:17) = nan;
        end
    end
end
boundaries = cellfun(@uint16,boundaries,'UniformOutput', 0);
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

%% Sinusoidal fit
function [amplitude, frequency, linear_length,arc_length] = sinusoidal(xy)

x = xy(:,1);
y = xy(:,2);

linear_length = sqrt((x(1)-x(end))^2 + (y(1)-y(end))^2);
arc_length = sum(sqrt((diff(x)).^2 + (diff(y)).^2));

p = polyfit(x,y,1);
R = [cos(atan(p(1))) -sin(atan(p(1))); sin(atan(p(1))) cos(atan(p(1)))]; %transformation matrix
t = [x,y]*R;
y_hat = t(:,2);
x_hat = t(:,1);

if mean(diff(x_hat))<0
    x_hat = flip(x_hat);
end

x_hat = x_hat-min(x_hat);
y_hat = y_hat-mean(y_hat);
Fs = 1/mean(diff(x_hat));
y_hat = interp1(x_hat, y_hat, min(x_hat):1/Fs:max(x_hat));

padding = 1000;
xdft = fft(y_hat, padding);
xdft = xdft(1:length(xdft)/2+1);
xdft = 1/length(y_hat).*xdft;
xdft(2:end-1) = 2*xdft(2:end-1);
freq = Fs/2*linspace(0,1,padding/2+1);
amplitude = max(abs(xdft));
frequency = freq(abs(xdft)==amplitude);

end

%% Arc length
function arc_length = find_arc_length(arc_length, boundary, extended)

p = arc_length(end-1:end,:);
if diff(p(:,1))<0
    b = diff(p(:,2))/diff(p(:,1));
    a = p(1,2)-b*p(1,1); %intersect
    x = flip(p(2)-extended:0.1:p(2));
    y = a + x*b;
elseif diff(p(:,1))>0
    b = diff(p(:,2))/diff(p(:,1));
    a = p(1,2)-b*p(1,1); %intersect
    x = p(2):0.1:p(2)+extended;
    y = a + x*b;
elseif diff(p(:,1))==0
    b = diff(p(:,1))/diff(p(:,2));
    a = p(1,1); %intersect
    if diff(p(:,2))>0
        y = p(4):0.1:p(4)+extended;
    else
        y = flip(p(4)-extended:0.1:p(4));
    end
    x = a + y*b;
end
arc_length = [arc_length;x(:),y(:)];

% arc_length(inpolygon(x',y',boundary(:,1),boundary(:,2)),:)
% row = row-extension;
% column = column-extension;
% arc_length = [arc_length;row',column'];

p = arc_length(1:2,:);
b = diff(p(:,2))/diff(p(:,1));
a = p(1,2)-b*p(1,1); %intersect
if diff(p(:,1))>0
    x = p(1)-extended:0.1:p(1);
    y = a + x*b;
elseif diff(p(:,1))<0
    x = flip(p(1):0.1:p(1)+extended);
    y = a + x*b;
elseif diff(p(:,1))==0
    b = diff(p(:,1))/diff(p(:,2));
    a = p(1,1); %intersect
    if diff(p(:,2))>0
        y = p(3)-extended:0.1:p(3);
    else
        y = flip(p(3):0.1:p(3)+extended);
    end
    x = a + y*b;
end
arc_length = [x',y';arc_length];

arc_length = arc_length(inpolygon(arc_length(:,1),arc_length(:,2),boundary(:,1),boundary(:,2)),:);

% p = flipud(arc_length(1:2,1:2));
% [row, column] = find_closest_edge(p, boundary, bw);
% arc_length = [row',column';arc_length];

end

%% Curvature
function [r,a,b, R2, df, RMSE] = radius_of_curvature(x,y)
%  given a load of points, with x,y coordinates, we can estimate the radius
%  of curvature by fitting a circle to them using least squares.
%    translate the points to the centre of mass coordinates
try
    mx = mean(x);
    my = mean(y);
    X = x - mx; Y = y - my;
    
    dx2 = mean(X.^2);
    dy2 = mean(Y.^2);
    
    %    Set up linear equation for derivative and solve
    RHS = (X.^2-dx2+Y.^2-dy2)/2;
    M = [X,Y];
    
    t = M\RHS;
    
    %    t is the centre of the circle [a0;b0]
    a0 = t(1); b0 = t(2);
    
    %   from which we can get the radius
    r = sqrt(dx2+dy2+a0^2+b0^2);
    
    %   return to given coordinate system
    a = a0 + mx;
    b = b0 + my;
    
    %% goodness of fit
    [curvefit,gof,~]=fit(x-a,(y-b).^2,@(beta,x) ((beta^2-x.^2)), 'Startpoint',r);
    
    if gof.adjrsquare<0||gof.adjrsquare>1, error('bad fit'),end
    
    r = curvefit.beta;
    R2 = gof.adjrsquare;
    df = gof.dfe;
    RMSE = gof.rmse;
    
catch
    r = nan;
    a = nan;
    b = nan;
    df = nan;
    R2 = nan;
    RMSE = nan;
end
end


