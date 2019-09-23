function [boundary_threshold, particle_threshold] =...
    threshold_level(V, background, maxV, minV, frames,...
    minimum_area, pixels_range, scaling, complement, noholes, ext)

% This function allows the user to interactively determine the boundary and
% particle thresholds of a given video.
%
% input variables:
%
%       V is the video matrix.
%
%       background is the average frame that is substracted from each frame
%           to get rid of non moving particles and noise.
%
%       maxV and minV are the limits of the de-backgrounded images
%           for histogram equalization.
%
%       frames is a vector with the list of frames to be used for
%           thresolding.
%
%       maxV and minV are the maximum and minimum intensity values of the
%           video after substracting the background and adding the mean
%           limits of the de-backgrounded images for histogram
%           equalization.
%
%       pixels_range is a 2X2 matrix with matrix indexes for the
%          portion of the video frames to be used in the analyses. Useful
%          to cut off information.
%
%       complement is a logical variable. If true, images are inverted
%           using imcomplement. Useful for light particles on dark field,
%           such as fluorescent images. 
%
%       noholes is a logic variable. If true, bwboundaries searches for
%           object (parent and child) boundaries. This can provide better
%           performance. If false, bwboundaries searches for both object
%           and hole boundaries. Default is false.
%
%
% output variables:
%
%       boundary_threshold (optional) is the minimum pixel intensity in a
%           grayscale (0-255) of the boundaries of a particle.
%
%       particle_thresholds (optional) is the minimum pixel intensity in a
%           grayscale (0-255) that a particle needs to have to be included
%           in the analyses.
%
%
% Called by particle_detection.m
%
% Calls bwboundaries_noholes, image_normalization
%
%       ex.: [boundary_threshold, particle_threshold] =...
%           threshold_level(V, background, maxV, minV, frames,...
%           minimum_area, pixels_range)
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
%  phase-contrast microscope

if strcmpi(ext,'.avi')
%     strcmp(class(V),'VideoReader') % video
    ff = frames(1);
    I = read(V,ff);
elseif size(V,3)>1 % stack of images
    ff = frames(1);
    I = V(:,:,ff);
elseif strcmpi(ext, '.tif') ||strcmpi(ext, '.tiff') 
%     isstruct(V) % multitiff
    ff = frames(1);
    I = imread(V(1).Filename,ff);
    if V(1).BitDepth==16
        if ~isempty(regexp(V(1).ImageDescription,'Hamamatsu', 'once'))
            I = uint8(bitshift(I,-4)); % in fact, when the V.Bitdepth =16 in a hamamatsu multitif, the bitdepth=12
        else
            I = uint8(I);
        end
    end
elseif strcmpi(ext, '.mat')
    ff = frames(1);
    I = imread([V.path_name, char(V.fns(ff))]);
else % single image
    I = V;
end

I = I(pixels_range(1,1):pixels_range(1,2),pixels_range(2,1):pixels_range(2,2));

if nargin<2 || isempty(background)
    background = zeros(size(I,1),size(I,2), 'uint8');
end

%% initialize outputs
if complement
%     background = imcomplement(background);
    I = imcomplement(I);
end
[pI, I] = image_normalization(I, background, minV, maxV, scaling);
% pI = imadjust(pI);
% m = uint8(mean(pI(:)));
m = uint8(graythresh(pI(:))*255);
boundary_threshold = m;
particle_threshold = boundary_threshold;
boundaries = {};

%% figure
f = figure('units','normalized','outerposition',[0 0 1 1]);
drawnow

set(f,'units', 'pixels')
p = get(f, 'position');
a = [p(3)-500 p(4)-350]./size(I);
if a(1)>a(2)
    hax(1) = axes('units', 'pixels');%,...
        set(hax(1),'position',[100 150 p(3)-300 (p(4)-250)*size(2)/size(1)]);
    Panel = uipanel('BackgroundColor','white',...
        'Units','pixels');%,...
        set(Panel,'Position',[p(3)-150 150 100 (p(4)-250)*size(2)/size(1)]);
else
    hax(1) = axes('units', 'pixels',...
        'position',[100 150 (p(3)-300)*size(1)/size(2) p(4)-200]);
    Panel = uipanel('BackgroundColor','white',...
        'Units','pixels',...
        'Position',[p(3)-150 150 100 p(4)-200]);
end
c = imshow(I,'parent',hax(1));
hold(hax(1), 'all')
h = plot(nan,nan,'b');
cd
slider(1) = uicontrol(Panel,'Style', 'slider',...
    'Min',1,'Max',255,'Value',boundary_threshold,...
    'units', 'normalized',...
    'Position', [0.40 0.05 0.25 0.95],...
    'BackgroundColor',[0.5 0.2 0.2],...
    'Callback', {@slide,true});

slider(2) = uicontrol(Panel,'Style', 'slider',...
    'Min',1,'Max',255,'Value',particle_threshold,...
    'SliderStep', [1/126 0.1],...
    'units', 'normalized',...
    'Position', [0.70 0.05 0.25 0.95],...
    'BackgroundColor',[0.2 0.2 0.5],...
    'Callback', {@slide,true});
chk = uicontrol(Panel,'Style','checkbox','String', 'Dark','units','normalized','Position', [0 0  1 0.05],'value',false,'Callback',@slide);
ticks = axes('Parent', Panel,'position',[0.40 0.07 0.01 0.907]);
set(ticks,'box','off','xtick',[],'xcolor',get(f,'Color'),'color','none','ylim',[1 255], 'TickDir','out')
if strcmp(class(V),'VideoReader') || isstruct(V)
    positions = get(hax(1),'position');
    slider(3) = uicontrol('Style', 'slider',...
        'Min', frames(1), 'Max', frames(end),...
        'SliderStep', [1/length(frames) 0.1],...
        'Value',ff,...
        'units', 'pixels',...
        'Position', [positions(1), 50 positions(3) 50],...
        'Callback', {@frame,hax(1)});
end

ok = uicontrol('style', 'pushbutton', ...
    'String', 'OK',...
    'units', 'pixels', ...
    'Position', [p(3)*2.25/20+p(4)*7/10*size(I,1)/size(I,2) p(4)/60 p(4)/30 p(4)/30],...
    'Callback', 'close');

slide(true)

uiwait(f)

%% change frame
    function frame(varargin)
        block_figure
        ff = round(get(slider(3), 'Value'));
        if strcmpi(ext,'.avi')
            %     strcmp(class(V),'VideoReader') % video
            I = read(V,ff);
        elseif size(V,3)>1 % stack of images
            ff = frames(1);
            I = V(:,:,ff);
            
        elseif strcmpi(ext, '.tif') ||strcmpi(ext, '.tiff')
            %     isstruct(V) % multitiff
            I = imread(V(1).Filename,ff);
            if V(1).BitDepth==16
                if ~isempty(regexp(V(1).ImageDescription,'Hamamatsu', 'once'))
                    I = uint8(bitshift(I,-4)); % in fact, when the V.Bitdepth =16 in a hamamatsu multitif,the bitdepth=12
                else
                    I = uint8(I);
                end
            end
        elseif strcmpi(ext, '.mat')
%             ff = frames(1);
            I = imread([V.path_name, char(V.fns(ff))]);
            
        end
        
        if complement
            I = imcomplement(I);
        end
        I = I(pixels_range(1,1):pixels_range(1,2),pixels_range(2,1):pixels_range(2,2));
        [pI, I] = image_normalization(I,background,minV,maxV, scaling);
%         pI = imadjust(pI);
        set(c,'Cdata', uint8(I)); % replots the image to the next frame
        find_boundaries
        unblock_figure
    end

    

%% slide function
    function slide(varargin)
        block_figure
        if get(chk,'value')
            boundary_threshold = get(slider(1),'Value');
            particle_threshold = get(slider(2),'Value');
            if particle_threshold>boundary_threshold
                particle_threshold = boundary_threshold;
                set(slider(2),'Value',particle_threshold)
            end
        else
            boundary_threshold = nan;
            particle_threshold = nan;
            boundaries = [];
        end
        find_boundaries
        unblock_figure
    end

    function find_boundaries(varargin)
        if get(chk,'value')
            boundaries = bwboundaries_noholes(pI, boundary_threshold,...
                particle_threshold, minimum_area, noholes); % identifies only particles without holes
        end
        replot
    end

    function replot(varargin)
        delete(h)
        h = [];
        for ii = 1:length(boundaries)
            h(ii) = plot(hax(1),boundaries{ii}(:,2),boundaries{ii}(:,1),'r');
        end
    end

    function block_figure(varargin)
        set(slider(:),'enable','off')
        set(chk,'enable','off')
        set(ok,'enable','off')
        pause(0.5)
    end

    function unblock_figure(varargin)
        if get(chk,'value')
            set(slider(1:2),'Enable','on');
        end
        set(slider(3),'Enable','on');
        set(chk,'enable','on')
        set(ok,'enable','on')
    end

end
