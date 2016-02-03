function [newdata, binVector] = binner(data,bincol,start,stop,binsize,threshold,bintype,fullvector)

% Bins a data array into uniform bin sizes using a user-defined bin type
% (e.g. @mean, @median, @var). Data array is assumed to be column oriented:
% rows = observations and columns = variables.
%
% Inputs:
%	data = data array to bin.
%	bincol = data column to use for setting the bins.
%	start = starting bin value.
%	stop = ending bin value.
%	binsize = bin size.
%   threshold is the minimum fraction of good values that a bin needs to have to
%       perform the operation. Otherwise, it will give NaN. For example, it the 
%       binsize is 2 and the threshold 0.5, those bins that have only one good
%       value will give NaN;
%	@bintype = function to use in calculating the bins.
%	fullvector. If 1 (conservative), it bins along the whole range defined bi
%       'start' and 'stop', producing nan values whenever there is no data to
%           bin.
%       If 0 (non-conservative), it bins only at the range of valid values
%           in data
%
% Outputs:
%	new = final, binned data array.
%
% Called by: moving_tracks
% 
% Examples:
%	new = bin([jd degC psu flr],1,jd(1),jd(end),1/60/24,@mean);
%	new = bin([db degC psu],1,2,100,2,@var);
%	new = bin(data,1,0,1,0.001,@median);
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
%  This file is part of trackbac, a collection of matlab scripts to
%  geometrically characterize and track swimming bacteria imaged by a
%  phase-contrast microscope.

if nargin==7
    fullvector=1;
end

% Removes rows with nans in the bincol
data(isnan(data(:,bincol)),:) = [];

% Set initial parameters
%%% bin vector and bin window
binVector = (start:binsize:stop)';
if isempty(binVector) || any(isnan(binVector))
    error('Start, Stop or binsize are not valid')
end
% if max(binVector)<max(data(:,bincol))
%     binVector(end+1) = binVector(end)+binsize;
% end

%determines the minimum number of valid values acceptable to perform the
%operation
fs = mode(diff(data(:,bincol)));
threshold = floor(binsize/fs*threshold);
if isinf(threshold)
    threshold=1;
end

if fullvector
    if max(data(:,bincol))<=binVector(end-1)
        data = [data; [binVector(end), nan(1,size(data,2)-1)]];
    end
else
    binVector = binVector(...
        binVector>=min(data(:,bincol)) &...
        binVector<=max(data(:,bincol)));
end

% removes values outside the range
data((data(:,bincol)<(binVector(1))|data(:,1)>(binVector(end))),:) = [];

% bin the data
[foo,b] = histc(data(:,bincol),binVector);
newdata = nan(size(binVector,1),size(data,2));
for ii = 1:size(data,2)
%     new(:,i) = accumarray(b, data(:,ii), [], bintype);
    newdata(:,ii) = accumarray(b(b>0),data(b>0,ii),size(binVector),bintype,nan);
end

%removes bins that have less than the minimum accepted threshold
for ii = 1:size(data,2)
    C = accumarray(b(b>0), data(b>0,ii),[],@(x) sum(~isnan(x)));
    newdata((C<threshold),ii) = nan;
end
newdata(:,bincol) = [];


%% alternatives
% this does the same but for large matrices is a little slower
% [~,b] = histc(data(:,bincol),binVector);
% b = [repmat(b(:),size(data,2),1) kron(1:size(data,2),ones(1,numel(b))).'];
% new = accumarray(b,data(:),[],bintype);
 
% [~,b] = histc(data(:,bincol),binVector);
% new = nan(nbins,size(data,2));
% for idx = 1:nbins
%     m = find(b==idx);
%     if isempty(m) == 1
%         continue
%     end
%     new(idx,:) = bintype(data(m,:),1);
% end
% %%% set the bin column equal to the binning vector
% new(:,bincol) = binVector;

% [~,b] = histc(data(:,bincol),binVector);
% for i = 1:nbins
%     C{i} = data(b==i,:);
% end %for
% new = cellfun(bintype,C,'UniformOutput',0);




