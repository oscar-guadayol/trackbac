function Runs_gui(fns)

% This gui allows to:
%   open a mat file or a list of mat files 
%   visualize individual tracks and their decomposition in runs and tumbles
%   call Runs and process all the tracks within the mat files
%   save the output from Runs.
%
% Input variables (optional):
%   fns: fullpath filename, a list of fullpath file names or a folder
%
% Calls: Runs
%
% Requires: SW_viscosity from Thermophysical properties of seawater toolbox
%        (http://web.mit.edu/seawater/)
%
% Ex.: Runs_gui('')
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

%% declare variables
pathname = pwd;
ext = [];
tracks = [];
ftracks = [];
if ~ exist('fns','var')
    fns = [];
end
FrameRate = nan;
scaling = nan;
win_length = 5;
methods_list = {'mean','median','max'};
method = 'mean';
tt = 1;
cell_width = 0.75*1e-6;      % mean width of cells in m
mean_cell_length = nan;
mean_cell_width = nan;
mean_cell_curvature = nan;
mean_cell_amplitude = nan;
mean_cell_frequency = nan;
time = nan;
speeds = nan;
angular_velocities = nan;
runs = nan;
tumbles = nan;
angle_threshold = [];
speed_threshold = [];
Geometries = [];
track = [];
track_lines = [];
moving = [];
percentage_moving = [];
run_lines = [];
tumble_lines = [];
Run_lengths = [];
Run_path_frequencies = [];
Run_speeds = [];
Run_complete = [];
Run_angle_deviations = [];
Run_mean_orientations = [];
Tumble_times = [];
Tumble_angles = [];
Tumble_angular_velocities = [];
Cell_area = [];
Cell_length = [];
Cell_width = [];
Cell_eccentricity = [];
Cell_equi_diam = [];
Cell_perimeter = [];
Cell_solidity = [];
Cell_radius = [];
Cell_radius_R2 = [];
Cell_radius_df = [];
Cell_radius_RMSE = [];
Cell_amplitude = [];
Cell_frequency = [];
StartTimes = [];
fig_tracks = [];
stop = false;
minimum_run_length = win_length/FrameRate;
% minimum_run_length = 0;

V = nan;
I = uint8(ones()*255);
h3d = figure; close(h3d)
fig_tracks = figure; close(fig_tracks)
dimensions_offset = 0.0; % offset bewteen the dimensions obtained with phase contrast microscopy and theoretical dimensions of standard beads
V = [];

%% figure
p = [1 27 1920 950];
f = figure('position', p);

% axis tracks
ax(1) = axes('position', [0.05 0.1 0.4 0.9]);

% plots first image of the first video
iptsetpref('ImshowAxesVisible','on');
p2 = imshow(I,'parent', ax(1));
xlabel(ax(1),'Distance (\mum)')
ylabel(ax(1),'Distance (\mum)')
set(ax(1),'ydir','normal');

hold(ax(1),'all')
c = plot(ax(1),nan, nan, '+y');
p = plot(ax(1), nan, nan, 'y-');
track_lines{1} = plot(ax(1), nan, nan, 'y-');

% axis angles
ax(2) = axes('position', [0.5 0.55 0.3 0.4]);
ylabel('Angular velocity (rad s^{-1})')
set(ax(2), 'yAxisLocation','right');
box(ax(2),'on');
hold(ax(2),'all');
angles_plot = plot(ax(2),nan,nan,'k');
hold(ax(2),'all')
timeline(2) = plot(ax(2),nan,'k');
angular_velocity_threshold = plot(ax(2),nan,nan,':k');

% axis speeds
ax(3) = axes('position', [0.5 0.1 0.3 0.4]);
ylabel('Speeds (\mum s^{-1})')
set(ax(3), 'yAxisLocation','right');
box(ax(3),'on');
hold(ax(3),'all');
speeds_plot = plot(ax(3), nan,nan,'k');
hold(ax(3),'all')
timeline(3) = plot(ax(3),nan,'k');
linear_velocity_threshold = plot(ax(3),nan,nan,':k');

%% panels
Panel.thresholds = uipanel('FontSize',12,  'title', 'Thresholds',...
    'Position',[0.86 0.95*12/19+0.03 0.13 0.95*7/19]);
Panel.Geometry = uipanel('FontSize',12,  'title', 'Geometry',...
    'Position',[0.86 0.95*7/19+0.03 0.13 0.95*5/19]);
Panel.control = uipanel('FontSize',12,  'title', 'Control',...
    'Position',[0.86 0.95*3/19+0.03 0.13 0.95*4/19]);
Panel.data= uipanel('FontSize',12,  'title', 'Data',...
    'Position',[0.86 0.03 0.13 0.95*3/19]);

%% uicontrols
rows = 8;
y_size = 1/(rows+(rows+1)/3);
y_positions = y_size/3*(1:rows)+y_size*(0:rows-1);
P(1) = uicontrol(Panel.thresholds,'Style','text','String','Window length (frames)',...
    'units','normalized','Position',[0.05 y_positions(8) 0.6 y_size], ...
    'TooltipString','Length of the smoothing filter window in seconds');
P(2) = uicontrol(Panel.thresholds,'Style','edit','String',win_length,...
    'units','normalized','Position',[0.7  y_positions(8) 0.25 y_size],...
    'Callback', {@stats, @thresholds}, ...
    'TooltipString','Length of the smoothing filter window in seconds');
P(3) = uicontrol(Panel.thresholds,'Style','text','String','Speed threshold (\mum/s)',...
    'units','normalized','Position',[0.05    y_positions(7)    0.6    y_size], ...
    'TooltipString','Minimum speed (\mum/s) of a run');
P(4) = uicontrol(Panel.thresholds,'Style','edit','String',speed_threshold,...
    'units','normalized','Position',[0.7    y_positions(7)    0.25    y_size],...
    'Callback', @thresholds, 'Enable', 'off',...
    'TooltipString','Minimum speed (\mum/s) of a run');
P(5) = uicontrol(Panel.thresholds,'Style','text','String','Angle threshold (rad/s)',...
    'units','normalized','Position',[0.05    y_positions(6)    0.6    y_size], ...
    'TooltipString','Minimum angle (rad/s) of a tumble');
P(6) = uicontrol(Panel.thresholds,'Style','edit','String',angle_threshold,...
    'units','normalized','Position',[0.7    y_positions(6)    0.25    y_size],...
    'Callback', @thresholds, 'Enable', 'off',...
    'TooltipString','Minimum angle (rad/s) of a tumble');
P(7) = uicontrol(Panel.thresholds,'Style','text','String','Minimum run length (s)',...
    'units','normalized','Position',[0.05   y_positions(5)    0.6    y_size], ...
    'TooltipString','Minimum length (\mum) of a run');
P(8) = uicontrol(Panel.thresholds,'Style','edit','String',minimum_run_length,...
    'units','normalized','Position',[0.7    y_positions(5)    0.25    y_size],...
    'Callback', @thresholds,...
    'TooltipString','Minimum length (\mum) of a run');
P(9) = uicontrol(Panel.thresholds,'Style','text','String','Smoothing mode',...
    'units','normalized','Position',[0.05   y_positions(4)    0.6    y_size], ...
    'TooltipString','Smoothing filter applied to tracks');
P(10) = uicontrol(Panel.thresholds,'Style','popup','String',methods_list ,...
    'units','normalized','Position',[0.7    y_positions(4)    0.25   y_size],...
    'Callback', @thresholds, ...
    'BusyAction', 'cancel', 'Interruptible', 'on');
P(24) =  uicontrol(Panel.thresholds,'Style','text','String','Dimensions offset (\mum)',...
    'units','normalized','Position',[0.05   y_positions(3)    0.6    y_size], ...
    'TooltipString','Smoothing filter applied to tracks');
P(25) = uicontrol(Panel.thresholds,'Style','edit','String',dimensions_offset,...
    'units','normalized','Position',[0.7    y_positions(3)    0.25    y_size],...
    'Callback', @thresholds,...
    'TooltipString','Dimension offset in a phase contrast microscope');
P(16) = uicontrol(Panel.thresholds,'Style','pushbutton','String','Histogram',...
    'units','normalized','Position',[0.05, y_positions(2), 0.9, y_size],...
    'Callback',@histogram3d, ...
    'BusyAction', 'cancel', 'Interruptible', 'on');
P(20) = uicontrol(Panel.thresholds,'Style','pushbutton','String','Export histogram',...
    'units','normalized','Position',[0.05, y_positions(1), 0.9, y_size],...
    'Callback',@export_histogram);

rows = 5;
y_size = 1/(rows+(rows+1)/3);
y_positions = y_size/3*(1:rows)+y_size*(0:rows-1);
P(26) = uicontrol(Panel.Geometry,'Style','text','String', 'Cell length (\mum)',...
    'units','normalized','Position',[0.05 y_positions(5) 0.6 y_size]);
P(27) = uicontrol(Panel.Geometry,'Style','edit','String', mean_cell_length,...
    'units','normalized','Position',[0.7  y_positions(5) 0.25 y_size], 'Enable', 'off');
P(28) = uicontrol(Panel.Geometry,'Style','text','String', 'Cell width (\mum)',...
    'units','normalized','Position',[0.05    y_positions(4)    0.6    y_size]);
P(29) = uicontrol(Panel.Geometry,'Style','edit','String', mean_cell_width,...
    'units','normalized','Position',[0.7    y_positions(4)    0.25    y_size],...
    'Enable', 'off');
P(30) = uicontrol(Panel.Geometry,'Style','text','String', 'Cell curvature',...
    'units','normalized','Position',[0.05    y_positions(3)    0.6    y_size]);
P(31) = uicontrol(Panel.Geometry,'Style','edit','String', mean_cell_curvature,...
    'units','normalized','Position',[0.7    y_positions(3)    0.25    y_size],...
    'Callback', @thresholds, 'Enable', 'off');
P(32) = uicontrol(Panel.Geometry,'Style','text','String', 'Helical amplitude',...
    'units','normalized','Position',[0.05   y_positions(2)    0.6    y_size]);
P(33) = uicontrol(Panel.Geometry,'Style','edit','String', mean_cell_amplitude,...
    'units','normalized','Position',[0.7    y_positions(2)    0.25    y_size],...
    'Enable', 'off');
P(34) = uicontrol(Panel.Geometry,'Style','text','String' ,'Helical frequency',...
    'units','normalized','Position',[0.05   y_positions(1)    0.6    y_size]);
P(35) = uicontrol(Panel.Geometry,'Style','edit','String', mean_cell_frequency ,...
    'units','normalized','Position',[0.7    y_positions(1)    0.25   y_size],...
    'Enable', 'off');

rows = 5;
y_size = 1/(rows+(rows+1)/3);
y_positions = y_size/3*(1:rows)+y_size*(0:rows-1);

Subpanel_tracks =  uipanel('Parent', Panel.control, 'title', 'Tracks',...
    'FontSize',12,  'units','normalized', 'Position',[0.05 y_positions(4) 0.9 y_size*7/3]);
P(11) = uicontrol(Subpanel_tracks, 'Style','Push','FontName','Blue Highway','String',char(171),...
    'FontSize',14,...
    'units','normalized','Position',[0.02 0.1 0.15 0.9], ...
    'Callback',  {@update_track,1}, ...
    'BusyAction', 'cancel', 'Interruptible', 'on');
P(12) = uicontrol(Subpanel_tracks, 'Style','Push','FontName','Blue Highway','String',char(60),...
    'FontSize',14,...
    'units','normalized','Position',[0.19 0.1 0.15 0.9], ...
    'Callback',  {@update_track,2}, ...
    'BusyAction', 'cancel', 'Interruptible', 'on');
P(13) = uicontrol(Subpanel_tracks,'Style','edit','String',tt,...
    'units','normalized','Position',[0.36 0.1 0.28 0.9],...
    'Callback', @thresholds, ...
    'BusyAction', 'cancel', 'Interruptible', 'on');
P(14) = uicontrol(Subpanel_tracks,'Style','Push','FontName','Blue Highway','String',char(62),...
    'FontSize',14,...
    'units','normalized','Position',[0.66 0.1 0.15 0.9], ...
    'Callback',  {@update_track,3}, ...
    'BusyAction', 'cancel', 'Interruptible', 'on');
P(15) = uicontrol(Subpanel_tracks, 'Style','Push','FontName','Blue Highway','String',char(187),...
    'FontSize',14,...
    'units','normalized','Position',[0.83 0.1 0.15 0.9], ...
    'Callback',  {@update_track,4}, ...
    'BusyAction', 'cancel', 'Interruptible', 'on');
P(17) = uicontrol(Panel.control,'Style','pushbutton','String','Replay',...
    'units','normalized','Position',[0.05, y_positions(3), 0.9, y_size],...
    'Callback',@thresholds, ...
    'BusyAction', 'cancel', 'Interruptible', 'on');
P(19) = uicontrol(Panel.control,'Style','pushbutton','String','Process all tracks',...
    'units','normalized','Position',[0.05, y_positions(2), 0.9, y_size],...
    'Callback',@process_all);
P(23) = uicontrol(Panel.control,'Style','pushbutton','String','Export figure',...
    'units','normalized','Position',[0.05, y_positions(1), 0.9, y_size],...
    'Callback',@export_figure);

rows = 3;
y_size = 1/(rows+(rows+1)/3);
y_positions = y_size/3*(1:rows)+y_size*(0:rows-1);
P(22) = uicontrol(Panel.data,'Style','pushbutton','String','Load data',...
    'units','normalized','Position',[0.05, y_positions(3), 0.9, y_size],...
    'Callback',{@open_mats,[]});
P(21) = uicontrol(Panel.data,'Style','pushbutton','String','Save data',...
    'units','normalized','Position',[0.05, y_positions(2), 0.9, y_size],...
    'Callback',@save_data);
P(18) = uicontrol(Panel.data,'Style','pushbutton','String','Close',...
    'units','normalized','Position',[0.05, y_positions(1), 0.9, y_size],...
    'Callback',@close_gui);

set(findobj('Type','uicontrol'), 'BusyAction','cancel', 'Interruptible','on');

if ~isempty(fns)
    open_mats(fns)
end

%% open data
    function open_mats(varargin)
        %% Data selection
        pause(0) % interrupts previous callbacks
        if isempty(varargin{3})
            fns = [];
        end
        if isempty(fns)
            [fns, pathname, ~] = uigetfile('multiselect','on',[pathname, filesep,  '*.mat']);
            fns = fullfile(pathname, fns);
        end % if isempty(fns)
        if ~iscell(fns)
            if isfolder(fns)
                [fns, pathname, ~] = uigetfile('multiselect','on',[fns, filesep,  '*.mat']);
                fns = fullfile(pathname, fns);
            end
        end % if isdir(fns)
        
        if ischar(fns)
            fns = {fns};
        end % if ischar(fns)
        
        if length(fns) == 1
            % looks for the variable fns inside the file, which indicates
            % that merge_tracks has been run on the dataset. If not, it
            % runs merge_tracks. Uses try/catch instead of whos because its
            % much faster.
            s = warning('error','MATLAB:load:variableNotFound');
            try
                m = matfile(char(fns));
                tracks = m.tracks;
                moving = m.moving;
                percentage_moving = m.percentage_moving;
                fns = m.fns;
                if ~iscell(fns)
                    [~,fns]=fileparts(char(m.fns));
                    fns = {[pathname fns]};
                end
            catch
                [tracks,moving, percentage_moving] =...
                    merge_tracks(fns, win_length);
            end
            warning(s)
        else
            [tracks,moving, percentage_moving] = merge_tracks(fns, win_length);
        end
        tracks = tracks(:,:,moving); % removes non moving cells
        
        clear Geometries
        for ff=1:length(fns)
            Geometries{ff} = load(fns{ff},'Boundaries', 'Centroids');
        end
                
        [pn, fn, ~] = fileparts(fns{tracks(1,end-1,1)});
        if ~strcmp(pn(end),filesep)
            pn = [pn, filesep];
        end
        data = load([pn, fn, '.mat'],'*Rate', 'scaling','fn');
        
        if ~isfield(data, 'FrameRate')
            FrameRate = data.SamplingRate;
        else
            FrameRate = data.FrameRate;
        end % if ~exist('FrameRate', 'var')
        scaling = data.scaling;

        [~, ~, ext] = fileparts(data.fn);
        if ~isempty(regexpi(ext, 'avi'))
            V = VideoReader([pn, fn, ext]);
            I = read(V,1);
        elseif ~isempty(regexpi(ext, '.tif*'))
            I = imread([pn, filesep, fn, '.tif'],1);
            if isa(I,'uint16')
                I = uint8(bitshift(I,-4)); % in fact, when the V.Bitdepth =16 in a hamamatsu multitif,the bitdepth=12
            end
            %% retrieves time information from multitiff
        elseif ~isempty(regexpi(ext, '.mat'))
            V = load([pn, fn, '_info']);
            V = V.V;
            I = imread([V.path_name char(V.fns(1))]);
        else
            I = ones(ceil(max(max((tracks(:,5,:))))*scaling),ceil(max(max((tracks(:,6,:))))*scaling))*255;
        end
        tt = 1;
        set(P(13),'String',tt);
        set(p2,'CData', I)
        axis(ax(1),'image')
        % formats xlabels to show \mums, not pixels
        xticks = (0:ceil(size(I,1)/scaling/10)*10/5:ceil(size(I,1)/scaling/10)*10);
        yticks = (0:ceil(size(I,2)/scaling/10)*10/5:ceil(size(I,2)/scaling/10)*10);
        xticklabels = reshape(sprintf('%4.0f', xticks),4,length(xticks))';
        yticklabels = reshape(sprintf('%4.0f', yticks),4,length(yticks))';
        set(ax(1),'xtick',xticks*scaling, 'xtickLabel',xticklabels)
        set(ax(1),'ytick',yticks*scaling, 'ytickLabel',yticklabels)
        
        minimum_run_length = win_length/FrameRate;
        set(P(8),'String', num2str(minimum_run_length))        
        
        stats
        thresholds
    end

%% update thresholds
    function thresholds(varargin)
        pause(0) % interrupts previous callbacks
        win_length = str2double(get(P(2),'String'));
        speed_threshold = str2double(get(P(4),'String'));
        angle_threshold = str2double(get(P(6),'String'));
        dimensions_offset = str2double(get(P(25),'String'));
        minimum_run_length = str2double(get(P(8),'String'));
        method = methods_list(get(P(10),'Value'));
        tt = str2double(get(P(13),'String'));
        stop = false;
        track_runs_and_tumbles
    end % if function update_track

%% update track number
    function update_track(varargin)
        pause(0) % interrupts previous callbacks
        if varargin{3}==3
            track_number = str2double(get(P(13),'String'))+1;
            if track_number<size(tracks,3)
                set(P(13),'String',num2str(track_number))
                thresholds
            end
        elseif varargin{3} == 2
            track_number = str2double(get(P(13),'String'))-1;
            if track_number>0
                set(P(13),'String',num2str(track_number))
                thresholds
            end
        elseif varargin{3} ==1
            track_number = 1;
            set(P(13),'String',track_number)
            thresholds
        elseif varargin{3} == 4
            track_number = size(tracks,3);
            set(P(13),'String',track_number)
            thresholds
        end %  if varargin{3}==3
    end % if function update_track

%% process all tracks
    function process_all(varargin)
        method = methods_list(get(P(10),'Value'));
        minimum_run_length = str2double(get(P(8),'String'));
        [Run_lengths, Run_speeds, Run_complete,...
            Run_angle_deviations, Run_mean_orientations, Run_path_frequencies, ...
            Tumble_times, Tumble_angles, Tumble_angular_velocities,...
            Cell_area, Cell_length, Cell_width, Cell_eccentricity, ...
            Cell_equi_diam, Cell_perimeter, Cell_solidity,...
            Cell_radius, Cell_radius_R2, Cell_radius_df, Cell_radius_RMSE,...
            Cell_amplitude, Cell_frequency,...
            StartTimes, fig_tracks] = ...
            Runs(tracks, cell_width, [], [], win_length,...
            minimum_run_length, method, FrameRate, dimensions_offset, 1, 0);
    end % function process_all

%% stats
    function stats(varargin)
        pause(0) % interrupts previous callbacks
        win_length = str2double(get(P(2),'String'));
        method = methods_list(get(P(10),'Value'));
        
        if win_length==1
            ftracks = tracks(:,5:6,:);
        elseif win_length>1
            ftracks = nan(size(tracks(:,5:6,:)));
            for ii = 1:size(tracks,3)
                ftracks(1:size(tracks(~isnan(tracks(:,5,ii)),5:6,ii),1),:,ii) =...
                    smooth_track(tracks(~isnan(tracks(:,5,ii)),5:6,ii), win_length, method);
            end % for ii
        end % if win_length==1
        x = diff(ftracks,1,1);
        
        speeds = sqrt(x(:,1,:).^2+x(:,2,:).^2)*FrameRate;
        speeds = speeds(:);

        Angles = atan2(x(:,2,:),x(:,1,:));
        angular_velocities = nan(size(Angles,1),size(Angles,3));
        angular_velocities(1:end-1,:) = abs(diff(unwrap(squeeze(Angles)),1,1));
        angular_velocities = angular_velocities(:);
        
    end % function stats

%% histogram
    function histogram3d(varargin)
        pause(0) % interrupts previous callbacks
        stats
        h3d = figure;
        hist3([speeds(~isnan(speeds)), angular_velocities(~isnan(speeds))],[100 100]);
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
        set(get(gca,'xLabel'),'String','Speeds (\mum/s)')
        set(get(gca,'yLabel'),'String','angular_velocities (Rad/s)')
        set(gca,'YDir','reverse')
    end % function histogram3d

%% export histogram
    function export_histogram(varargin)
        pause(0) % interrupts previous callbacks
        if ~isvalid(h3d)
            return
        end % if isempty(fig_tracks)
        load plot_style
        plot_style.Color= 'rgb';
        plot_style.Width = 21.40;
        plot_style.Height = 16.84;
        plot_style.Resolution = 72;
        plot_style.FontName = 'Helvetica';
        plot_style.FixedFontSize = '16';
        plot_style.Format = 'eps';
        
        %default folder and filename
        fn = double(char(fns));
        fn = [char(fn(1,all(bsxfun(@eq, fn, fn(1,:))))),'.jpeg'];
        [fn,pn] = uiputfile('*.jpeg;*.eps','Export figure', fn);
        [~, ~, ext] = fileparts(fn);
        plot_style.Format = ext(2:end);
        figname = fullfile(pn, fn);
        
        % export figure
        hgexport(h3d,figname,plot_style)
    end % function export_figure

%% export figure
    function export_figure(varargin)
        pause(0) % interrupts previous callbacks
        if ~isvalid(fig_tracks)
            return
        end % if isempty(fig_tracks)
        load plot_style
        plot_style.Color= 'rgb';
        plot_style.Width = 21.40;
        plot_style.Height = 16.84;
        plot_style.Resolution = 72;
        plot_style.FontName = 'Helvetica';
        plot_style.FixedFontSize = '16';
        plot_style.Format = 'eps';
        
        %default folder and filename
        fn = double(char(fns));
        fn = [char(fn(1,all(bsxfun(@eq, fn, fn(1,:))))),'.jpeg'];
        [fn,pn] = uiputfile('*.jpeg;*.eps','Export figure', fn);
        [~, ~, ext] = fileparts(fn);
        plot_style.Format = ext(2:end);
        figname = fullfile(pn, fn);
        
        % export figure
        hgexport(fig_tracks,figname,plot_style)
    end % function export_figure

%% save data
    function save_data(varargin)
        % default folder and filename
        pause(0) % interrupts previous callbacks
        fn = double(char(fns));
        fn = char(fn(1,all(bsxfun(@eq, fn, fn(1,:)))));
        [pn,fn] = fileparts(fn);
        [fn,pn] = uiputfile('*.mat;*.csv','Save data', [pn, filesep, fn, '_runs.mat']);
        fn = fullfile(pn, fn);
        [~, ~, ext] = fileparts(fn);
        if strcmp(ext, '.mat')
            save(fn, 'Run_lengths', 'Run_speeds','Run_complete',...
                'Run_angle_deviations','Run_mean_orientations', 'Run_path_frequencies',...
                'Tumble_times', 'Tumble_angles', 'Tumble_angular_velocities',...
                'Cell_area', 'Cell_length', 'Cell_width',...
                'Cell_eccentricity', 'Cell_equi_diam','Cell_perimeter', 'Cell_solidity',...
                'Cell_radius', 'Cell_radius_R2', 'Cell_radius_df', 'Cell_radius_RMSE',...
                'Cell_amplitude', 'Cell_frequency',...
                'win_length',...
                'minimum_run_length','StartTimes', 'tracks', 'moving', 'percentage_moving', '-v7.3');
        elseif strcmp(ext, '.csv')
             Headers = {'Track_number',',',...
                 'Run_length (s)',',',...
                 'Mean_run_speed (micron/s^',',',...
                 'Std_Run_speed (micron/s)',',',...
                 'Min_run_speed (micron/s)',',',...
                 'Max_run_speed (micron/s)',',',...
                 'Tumble_times (s)',',',...
                 'Tumble_angles (radians)',',',...
                 'Mean_cell_area (micron^2)',',',...
                 'Std_cell_area (micron^2)',',',...
                 'Mean_Cell_length (micron)',',',...
                 'Std_Cell_length (micron)',',',...
                 'Mean_Cell_width (micron)',',',...
                 'Std_Cell_width (micron)',',',...
                 'Mean_Cell_eccentricity',',',...
                 'Std_Cell_eccentricity',',',...
                 'Mean_Cell_equi_diam',',',...
                 'Std_Cell_equi_diam',',',...
                 'Mean_Cell_perimeter (micron)',',',...
                 'Std_Cell_perimeter (micron)',',',...
                 'Mean_Cell_solidity',',',...
                 'Std_Cell_solidity',',',...
                 'Mean_Cell_radius (micron)',',',...
                 'Std_Cell_radius (micron)',',',...
                 'Mean_Cell_radius_R2',',',...
                 'Std_Cell_radius_R2',',',...
                 'Mean_Cell_radius_df',',',...
                 'Std_Cell_radius_df',',',...
                 'Mean_Cell_radius_RMSE',',',...
                 'Std_Cell_radius_RMSE',',',...
                 'Mean_Cell_amplitude (micron)',',',...
                 'Std_Cell_amplitude (micron)',','...
                 'Mean_Cell_frequency (micron^-1)',',',...
                 'Std_Cell_frequency (micron^-1)',','};
             Headers = cell2mat(Headers);
            D = [Run_lengths(:,1),Run_lengths(:,2),...
                Run_speeds(:,2), Run_speeds(:,3), Run_speeds(:,4), Run_speeds(:,5),...
                Tumble_times, Tumble_angles,...
                Cell_area(:,1),Cell_area(:,2),...
                Cell_length(:,1),Cell_length(:,2),...
                Cell_width(:,1),Cell_width(:,2),...
                Cell_eccentricity(:,1),Cell_eccentricity(:,2),...
                Cell_equi_diam(:,1), Cell_equi_diam(:,2),...
                Cell_perimeter(:,1),Cell_perimeter(:,2),...
                Cell_solidity(:,1),Cell_solidity(:,2),...
                Cell_radius(:,1),Cell_radius(:,2),...
                Cell_radius_R2(:,1),Cell_radius_R2(:,2),...
                Cell_radius_df(:,1),Cell_radius_df(:,2),...
                Cell_radius_RMSE(:,1),Cell_radius_RMSE(:,2),...
                Cell_amplitude(:,1),Cell_amplitude(:,2),...
                Cell_frequency(:,1),Cell_frequency(:,2)];            
            fileID = fopen(fn,'w');
            fprintf(fileID,'%s\n',Headers);
            fclose(fileID);
            dlmwrite(fn, D, 'precision', '%10.2f', '-append');
            
        end % if strcmp(ext, '.mat')
    end % function save_data

%% runs and tumbles
    function track_runs_and_tumbles(varargin)
        pause(0) % interrupts previous callbacks
        method = methods_list(get(P(10),'Value'));
        minimum_run_length = str2double(get(P(8),'String'));
        track = tracks(:,:,tt);
        a = (nanmedian(track(:,15))-dimensions_offset)*1e-6;        
        b = (nanmedian(track(:,17))-dimensions_offset)*1e-6;
        mean_cell_length = [nanmedian(track(:,15))-dimensions_offset nanstd(track(:,15))-dimensions_offset];
        mean_cell_width = [nanmedian(track(:,17))-dimensions_offset nanstd(track(:,17))-dimensions_offset];
        mean_cell_curvature = [nanmean(track(:,18)) nanstd(track(:,18))];
        mean_cell_amplitude = [nanmean(track(:,22)) nanstd(track(:,22))];
        mean_cell_frequency = [nanmean(track(:,23)) nanstd(track(:,23))];
        
        P(27).String = sprintf(['%0.2f' ' ' char(177) ' ' '%0.2f'], mean_cell_length);
        P(29).String = sprintf(['%0.2f' ' ' char(177) ' ' '%0.2f'], mean_cell_width);
        P(31).String = sprintf(['%0.2f' ' ' char(177) ' ' '%0.2f'], mean_cell_curvature);
        P(33).String = sprintf(['%0.2f' ' ' char(177) ' ' '%0.2f'], mean_cell_amplitude);
        P(35).String = sprintf(['%0.2f' ' ' char(177) ' ' '%0.2f'], mean_cell_frequency);        

        time = (0:1/FrameRate:(size(track,1)-1)/FrameRate)';
        time(isnan(track(:,5)),:) = [];
        track(isnan(track(:,5)),:) = [];
        
        %% Brownian thresholds
        [~,~, Dt, Dr]  = theoretical_friction_coefficients(a/2,cell_width/2,cell_width/2,30); 
        p_value = 1-99.7300204/100;
        time_step = win_length/FrameRate;
        rms_angle = (2*Dr(2)*time_step).^(1/2); % standard deviation of the normal distribution of angles after a time=time_step
        angle = abs(norminv(p_value/2,0,rms_angle));
        rms_displacement = (2*Dt(1)*time_step).^(1/2); % standard deviation of the normal distribution of displacements after a time=time_step
        displacement = norminv(1-p_value/2,0,rms_displacement);
        if isempty(angle_threshold) || isnan(angle_threshold)
            angle_threshold = angle/time_step; % apparent angular velocity threshold in radians per second
        end % if isempty(angle_threshold)
        if isempty(speed_threshold) || isnan(speed_threshold)
            speed_threshold = displacement*1e6/time_step; % apparent speed threshold in µmeters per second during time step
        end % if isempty(speed_threshold)
        
        [~, ~, ~, ~, ~, ~, ~, ~, ~, speeds, angular_velocities, runs, tumbles,...
            time, ~, angle_threshold, speed_threshold] =...
            runs_and_tumbles(track, FrameRate, angle_threshold, speed_threshold,...
            minimum_run_length, a, b, method, win_length);
        
        P(4).String = num2str(speed_threshold);
        P(6).String = num2str(angle_threshold);
        track_plot % calls ploting
    end % function runs_and_tumbles(varargin)

%% closes gui
    function close_gui(varargin)
        pause(0) % interrupts previous callbacks
        close(f)
    end

%% Track_run
    function track_plot(varargin)
        pause(0) % interrupts previous callbacks
        %% loads right video
        [pn, fn, ~] = fileparts(fns{tracks(1,end-1,tt)});
        if ~isempty(regexpi(ext, 'avi'))
            V = VideoReader([pn, filesep, fn, ext]);
            I = read(V, track(1,2));
        elseif ~isempty(regexpi(ext, '.tif'))
            I = imread([pn, filesep, fn, '.tif'],track(1,2));
        elseif ~isempty(regexpi(ext, '.mat'))
            I = imread([V.path_name char(V.fns(track(1,2)))]);
        end
        
        %%
        set(p2,'CData', I)
        
        set(angles_plot,'xdata',time,'ydata',angular_velocities)
        set(timeline(2),'xdata',nan,'ydata',nan)
        angular_velocity_threshold.XData = time; 
        angular_velocity_threshold.YData = angle_threshold*ones(size(time));
        
        linear_velocity_threshold.XData = time; 
        linear_velocity_threshold.YData = speed_threshold*ones(size(time));
        
        for ii =1:length(tumble_lines)
            set(tumble_lines(ii),'xData',nan,'yData',nan)
        end
        for ii =1:length(run_lines)
            set(run_lines(ii),'xData',nan,'yData',nan)
        end
                
        for  ii = 1 : size(tumbles,1)
            tumble_lines(ii) = plot(ax(2),time(tumbles(ii,1)-1:tumbles(ii,2)+1),angular_velocities(tumbles(ii,1)-1:tumbles(ii,2)+1),'+-b');
        end
        
        for ii=1:size(runs,1)
            run_lines(ii) = plot(ax(2),time(runs(ii,1):runs(ii,2)),angular_velocities(runs(ii,1):runs(ii,2)),'+-r');
        end
        
        set(speeds_plot,'xdata',time,'ydata',speeds)
        set(timeline(3),'xdata',nan,'ydata',nan)
        
        for ii = 1:size(tumbles,1)
            tumble_lines(ii + size(tumbles,1)) = plot(ax(3),time(tumbles(ii,1)-1:tumbles(ii,2)+1),speeds(tumbles(ii,1)-1:tumbles(ii,2)+1),'+-b');
        end
        
        for ii=1:size(runs,1)
            run_lines(ii + size(runs,1)) = plot(ax(3),time(runs(ii,1):runs(ii,2)),speeds(runs(ii,1):runs(ii,2)),'+-r');
        end
        
        
        for ll = 1:length(track_lines)
            set(track_lines{ll},'xdata', nan, 'ydata', nan)
        end
        ll=1;
        
        for ff = 1:size(track,1) %frames
            if stop
                set(timeline(2),'xdata',nan,'ydata',nan)
                set(timeline(3),'xdata',nan,'ydata',nan)
                break
            end
            %             try
            timerVal(1) = tic;
            
            if ~isempty(regexpi(ext, 'avi'))
                V = VideoReader([pn, filesep, fn, ext]);
                I = read(V, track(ff,2));
            elseif ~isempty(regexpi(ext, '.tif'))
                I = imread([pn, filesep, fn, '.tif'], track(ff,2));
                if isa(I,'uint16')
                    I = uint8(bitshift(I,-4)); % in fact, when the V.Bitdepth =16 in a hamamatsu multitif,the bitdepth=12
                end
            elseif ~isempty(regexpi(ext, '.mat'))
                I = imread([V.path_name char(V.fns(track(ff,2)))]);                
            end
                    
            set(p2,'Cdata',I)
            set(p, 'Xdata', nan, 'YData', nan)
            set(c, 'Xdata', nan, 'YData', nan)
            
            if ~isempty(Geometries{track(1,end-1)}.Boundaries)            
                x = Geometries{track(1,end-1)}.Boundaries{track(ff,2)}{track(ff,3)}(:,2);
                y = Geometries{track(1,end-1)}.Boundaries{track(ff,2)}{track(ff,3)}(:,1);
                set(p, 'Xdata', x, 'YData', y)
                set(c, 'Xdata', Geometries{track(1,end-1)}.Centroids{track(ff,2)}(track(ff,3),1),...
                       'YData', Geometries{track(1,end-1)}.Centroids{track(ff,2)}(track(ff,3),2))
            end
             set(track_lines{1}, 'Xdata', track(1:ff,5)*scaling,...
                    'YData', track(1:ff,6)*scaling)
            

            if ~isempty(tumbles) && ismember(ff,tumbles(:,1))
                ff0 = ff;
                ll = ll+1;
                color = 'b';
                track_lines{ll} = plot(ax(1), nan, nan, 'color',color);
            end
            if ~isempty(runs) && ismember(ff,runs(:,1))
                ff0 = ff;
                ll = ll+1;
                color = 'r';
                track_lines{ll} = plot(ax(1), nan, nan, 'color',color);
            end
            
            if (~isempty(runs) || ~isempty(tumbles)) && ff>max([runs(:,2); tumbles(:,2)])
                ll = 1;
            end
            
            if ll>1
                if ff0 == 1, ff0 = 2; end
                set(track_lines{ll}, 'Xdata', track(ff0-1:ff,5)*scaling,...
                    'YData', track(ff0-1:ff,6)*scaling)
            end
            
            set(timeline(2),'xdata',[time(ff) time(ff)],'ydata',get(ax(2), 'ylim' ))
            set(timeline(3),'xdata',[time(ff) time(ff)],'ydata',get(ax(3), 'ylim' ))
            if 1/FrameRate>toc(timerVal)
                pause(1/FrameRate-toc(timerVal))
            end
            drawnow
        end
        stop = true;
    end
end