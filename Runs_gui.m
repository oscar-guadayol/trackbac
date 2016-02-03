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
tracks = nan(1,6,1);
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
% b = 1.3374*1e-6; % mean width of cells in m
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
Run_speeds = [];
Run_angle_deviation = [];
Run_cos_deviation = [];
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
StartTimes = [];
fig_tracks = [];
stop = false;
minimum_run_length = win_length/FrameRate;
V = nan;
I = uint8(ones()*255);
h3d = figure; close(h3d)
fig_tracks = figure; close(fig_tracks)

%% figure
p = [1 27 1920 950];
f = figure('position', p);

% axis tracks
ax(1) = axes('position', [0.05 0.1 0.4 0.9]);

% plots first image of the first video
iptsetpref('ImshowAxesVisible','on');
p2 = imshow(I,'parent', ax(1));
xlabel(ax(1),'Distance (µm)')
ylabel(ax(1),'Distance (µm)')
set(ax(1),'ydir','normal');

hold(ax(1),'all')
c = plot(ax(1),nan, nan, '+y');
p = plot(ax(1), nan, nan, 'y-');
track_lines{1} = plot(ax(1), nan, nan, 'y-');

% axis angles
ax(2) = axes('position', [0.5 0.55 0.3 0.4]);
ylabel('\deltaAngle / \deltat (rad s^{-1})')
set(ax(2), 'yAxisLocation','right');
box(ax(2),'on');
hold(ax(2),'all');
angles_plot = plot(ax(2),nan,nan,'k');
hold(ax(2),'all')
timeline(2) = plot(ax(2),nan,'r');

% axis speeds
ax(3) = axes('position', [0.5 0.1 0.3 0.4]);
ylabel('Speeds (µm s^{-1})')
set(ax(3), 'yAxisLocation','right');
box(ax(3),'on');
hold(ax(3),'all');
speeds_plot = plot(ax(3), nan,nan,'k');
hold(ax(3),'all')
timeline(3) = plot(ax(3),nan,'r');

%% panels
Panel.thresholds = uipanel('FontSize',12,  'title', 'Thresholds',...
    'Position',[0.86 0.9*7/14+0.05 0.13 0.9*7/14]);
Panel.control = uipanel('FontSize',12,  'title', 'Control',...
    'Position',[0.86 0.9*3/14+0.05 0.13 0.9*4/14]);
Panel.data= uipanel('FontSize',12,  'title', 'Data',...
    'Position',[0.86 0.05 0.13 0.9*3/14]);

%% uicontrols
rows = 7;
y_size = 1/(rows+(rows+1)/3);
y_positions = y_size/3*(1:rows)+y_size*(0:rows-1);
P(1) = uicontrol(Panel.thresholds,'Style','text','String','Window length (frames)',...
    'units','normalized','Position',[0.05 y_positions(7) 0.6 y_size], ...
    'TooltipString','Length of the smoothing filter window in seconds');
P(2) = uicontrol(Panel.thresholds,'Style','edit','String',win_length,...
    'units','normalized','Position',[0.7  y_positions(7) 0.25 y_size],...
    'Callback', {@stats, @thresholds}, ...
    'TooltipString','Length of the smoothing filter window in seconds');
P(3) = uicontrol(Panel.thresholds,'Style','text','String','Speed threshold (µm/s)',...
    'units','normalized','Position',[0.05    y_positions(6)    0.6    y_size], ...
    'TooltipString','Minimum speed (µm/s) of a run');
P(4) = uicontrol(Panel.thresholds,'Style','edit','String',speed_threshold,...
    'units','normalized','Position',[0.7    y_positions(6)    0.25    y_size],...
    'Callback', @thresholds, 'Enable', 'off',...
    'TooltipString','Minimum speed (µm/s) of a run');
P(5) = uicontrol(Panel.thresholds,'Style','text','String','Angle threshold (rad/s)',...
    'units','normalized','Position',[0.05    y_positions(5)    0.6    y_size], ...
    'TooltipString','Minimum angle (rad/s) of a tumble');
P(6) = uicontrol(Panel.thresholds,'Style','edit','String',angle_threshold,...
    'units','normalized','Position',[0.7    y_positions(5)    0.25    y_size],...
    'Callback', @thresholds, 'Enable', 'off',...
    'TooltipString','Minimum angle (rad/s) of a tumble');
P(7) = uicontrol(Panel.thresholds,'Style','text','String','Minimum run length (s)',...
    'units','normalized','Position',[0.05   y_positions(4)    0.6    y_size], ...
    'TooltipString','Minimum length (µm) of a run');
P(8) = uicontrol(Panel.thresholds,'Style','edit','String',minimum_run_length,...
    'units','normalized','Position',[0.7    y_positions(4)    0.25    y_size],...
    'Callback', @thresholds, 'Enable', 'off',...
    'TooltipString','Minimum length (µm) of a run');
P(9) = uicontrol(Panel.thresholds,'Style','text','String','Smoothing mode',...
    'units','normalized','Position',[0.05   y_positions(3)    0.6    y_size], ...
    'TooltipString','Smoothing filter applied to tracks');
P(10) = uicontrol(Panel.thresholds,'Style','popup','String',methods_list ,...
    'units','normalized','Position',[0.7    y_positions(3)    0.25   y_size],...
    'Callback', @thresholds, ...
    'BusyAction', 'cancel', 'Interruptible', 'on');
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
        if isempty(varargin{3})
            fns = [];
        end
        if isempty(fns)
            [fns, pathname, ~] = uigetfile('multiselect','on',[pathname, filesep,  '*.mat']);
            fns = fullfile(pathname, fns);
        end % if isempty(fns)
        if ~iscell(fns)
            if isdir(fns)
                [fns, pathname, ~] = uigetfile('multiselect','on',[fns, filesep,  '*.mat']);
                fns = fullfile(pathname, fns);
            end
        end % if isdir(fns)
        
        if ischar(fns)
            fns = {fns};
        end % if ischar(fns)
        
        if length(fns) == 1
            ismerged = whos('-file',char(fns));
            if ismember('fns', {ismerged.name})
                tracks =  load(char(fns), 'tracks');
                tracks = tracks.tracks;
                moving =  load(char(fns), 'moving');
                moving = moving.moving;
                percentage_moving =  load(char(fns), 'percentage_moving');
                percentage_moving = percentage_moving.percentage_moving;
                fns =  load(char(fns), 'fns');
                fns = fns.fns;
            else
                [tracks,moving, percentage_moving] = merge_tracks(fns);
            end
        else
            [tracks,moving, percentage_moving] = merge_tracks(fns);
        end
        data = load(fns{1},'*Rate', 'scaling');
        if ~isfield(data, 'FrameRate')
            FrameRate = data.SamplingRate;
        else
            FrameRate = data.FrameRate;
        end % if ~exist('FrameRate', 'var')
        scaling = 2.75;
        clear data
        
        load(fns{tracks(1,end-1,1)}, 'fn')
        ext = fn(end-2:end);
        [pn, fn, ~] = fileparts(fns{tracks(1,end-1,1)});
        if strcmpi(ext,'avi')
            V = VideoReader([pn, filesep, fn, '.avi']);
            I = read(V,1);
        elseif strcmpi(ext, 'tif')
            V = imfinfo([pn, filesep, fn, '.tif']);
            I = imread(V(1).Filename,1);
            %% retrieves time information from multitiff
        end
        tt = 1;
        set(P(13),'String',tt);
        set(p2,'CData', I)
        axis(ax(1),'image')
        % formats xlabels to show µms, not pixels
        xticks = (0:ceil(size(I,1)/scaling/10)*10/5:ceil(size(I,1)/scaling/10)*10);
        yticks = (0:ceil(size(I,2)/scaling/10)*10/5:ceil(size(I,2)/scaling/10)*10);
        xticklabels = reshape(sprintf('%3.0f', xticks),3,length(xticks))';
        yticklabels = reshape(sprintf('%3.0f', yticks),3,length(yticks))';
        set(ax(1),'xtick',xticks*scaling, 'xtickLabel',xticklabels)
        set(ax(1),'ytick',yticks*scaling, 'ytickLabel',yticklabels)
        stats
        thresholds
    end

%% update thresholds
    function thresholds(varargin)
        win_length = str2double(get(P(2),'String'));
        speed_threshold = str2double(get(P(4),'String'));
        angle_threshold = str2double(get(P(6),'String'));
        minimum_run_length = win_length/FrameRate;
        set(P(8),'String', num2str(minimum_run_length))
        method = methods_list(get(P(10),'Value'));
        tt = str2double(get(P(13),'String'));
        stop = false;
        track_runs_and_tumbles
    end % if function update_track

%% update track number
    function update_track(varargin)
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
        [Run_lengths, Run_speeds, Run_angle_deviation, Run_cos_deviation,...
            Tumble_times, Tumble_angles, Tumble_angular_velocities,...
            Cell_area, Cell_length, Cell_width, Cell_eccentricity, ...
            Cell_equi_diam, Cell_perimeter, Cell_solidity,...
            StartTimes, fig_tracks] = ...
            Runs(tracks(:,:,:), fns, [], [], win_length,...
            minimum_run_length, method, FrameRate, 0, 1, 0);
    end % function process_all

%% stats
    function stats(varargin)
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
        speeds = speeds(:,:,moving);
        speeds = speeds(:);

        Angles = atan2(x(:,2,:),x(:,1,:));
        Angles(:,:,moving);
        angular_velocities = nan(size(Angles,1),size(Angles,3));
        angular_velocities(1:end-1,:) = abs(diff(unwrap(squeeze(Angles)),1,1));
        angular_velocities = angular_velocities(:);
        
    end % function stats

%% histogram
    function histogram3d(varargin)
        stats
        h3d = figure;
        hist3([speeds(~isnan(speeds)), angular_velocities(~isnan(speeds))],[100 100]);
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
        set(get(gca,'xLabel'),'String','Speeds (µm/s)')
        set(get(gca,'yLabel'),'String','angular_velocities (Rad/s)')
        set(gca,'YDir','reverse')
    end % function histogram3d

%% export histogram
    function export_histogram(varargin)
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
        fn = double(char(fns));
        fn = char(fn(1,all(bsxfun(@eq, fn, fn(1,:)))));
        [pn,fn] = fileparts(fn);
        [fn,pn] = uiputfile('*.mat;*.csv','Save data', [pn, filesep, fn, '_runs.mat']);
        fn = fullfile(pn, fn);
        [~, ~, ext] = fileparts(fn);
        if strcmp(ext, '.mat')
            save(fn, 'Run_lengths', 'Run_speeds', 'Run_angle_deviation',...
                'Tumble_times', 'Tumble_angles', 'Tumble_angular_velocities',...
                'Cell_area', 'Cell_length', 'Cell_width',...
                'Cell_eccentricity', 'Cell_equi_diam','Cell_perimeter', 'Cell_solidity',...
                'angle_threshold','speed_threshold', 'win_length',...
                'minimum_run_length','StartTimes', 'tracks', 'moving', 'percentage_moving', '-v7.3')
        elseif strcmp(ext, '.csv')
            rnames = {'Track','Run_length',...
                'Mean_run_speed','Std_Run_speed','Min_run_speed','Max_run_speed',...
                'Tumble_times', 'Tumble_angles',...
                'Mean_cell_area', 'Std_cell_area',...
                'Mean_Cell_length','Cell_length',...
                'Mean_Cell_width','Std_Cell_width',...
                'Mean_Cell_eccentricity','Std_Cell_eccentricity',...
                'Mean_Cell_equi_diam','Std_Cell_equi_diam',...
                'Mean_Cell_perimeter','Std_Cell_perimeter',...
                'Mean_Cell_solidity','Std_Cell_solidity'};
            D = dataset(Run_lengths(:,1),Run_lengths(:,2), Run_speeds(:,1),...
                Run_speeds(:,2), Run_speeds(:,3),...
                Run_speeds(:,4), Tumble_times, Tumble_angles, Cell_area(:,1),Cell_area(:,2),Cell_length(:,1),Cell_length(:,2),...
                Cell_width(:,1),Cell_width(:,2),Cell_eccentricity(:,1),Cell_eccentricity(:,2),Cell_equi_diam(:,1),...
                Cell_equi_diam(:,2),Cell_perimeter(:,1),Cell_perimeter(:,2),Cell_solidity(:,1),Cell_solidity(:,2), ...
                'Varnames', rnames);
            export(D,'file',fn)
        end % if strcmp(ext, '.mat')
    end % function save_data

%% runs and tumbles
    function track_runs_and_tumbles(varargin)
        method = methods_list(get(P(10),'Value'));
        minimum_run_length = str2double(get(P(8),'String'));
        track = tracks(:,:,tt);
        Geometries = load(fns{tracks(1,end-1,tt)},'Boundaries', 'Centroids');
        dimensions_offset = 0; % offset bewteen the dimensions obtained with phase contrast and those from electon microscope
        a = (nanmedian(track(:,15))-dimensions_offset)*1e-6;
%         b = (nanmedian(track(:,17))-dimensions_offset)*1e-6;
        b = 0.7*1e-6;
        
        track = tracks(:,:,tt);
        time = (0:1/FrameRate:(size(track,1)-1)/FrameRate)';
        time(isnan(track(:,5)),:) = [];
        track(isnan(track(:,6)),:) = [];
        Geometries = load(fns{tracks(1,end-1,tt)},'Boundaries', 'Centroids');
        [~, ~, ~, ~, ~, ~, ~, speeds, angular_velocities, runs, tumbles, time] =...
            runs_and_tumbles(track, FrameRate, angle_threshold, speed_threshold, minimum_run_length, a, b, method, win_length, '', 0, 0);                
        track_plot % calls ploting        
    end % function runs_and_tumbles(varargin)

%% closes gui
    function close_gui(varargin)
        close(f)
    end

%% Track_run
    function track_plot(varargin)
        
        %% loads right video
        [pn, fn, ~] = fileparts(fns{tracks(1,end-1,tt)});
        if strcmpi(ext,'avi')
            V = VideoReader([pn, filesep, fn, '.avi']);
            I = read(V,track(1,2));
        elseif strcmpi(ext, 'tif')
            V = imfinfo([pn, filesep, fn, '.tif']);
            I = imread(V(1).Filename, track(1,2));
        end
        
        %%
        set(p2,'CData', I)
        
        set(angles_plot,'xdata',time,'ydata',angular_velocities)
        set(timeline(2),'xdata',nan,'ydata',nan)
        
        for ii =1:length(tumble_lines)
            set(tumble_lines(ii),'xData',nan,'yData',nan)
        end
        for ii =1:length(run_lines)
            set(run_lines(ii),'xData',nan,'yData',nan)
        end
        
        for ii=1:size(runs,1)
            run_lines(ii) = plot(ax(2),time(runs(ii,1):runs(ii,2)),angular_velocities(runs(ii,1):runs(ii,2)),'r');
        end
        
        for  ii = 1 : size(tumbles,1)
            tumble_lines(ii) = plot(ax(2),time(tumbles(ii,1)-1:tumbles(ii,2)+1),angular_velocities(tumbles(ii,1)-1:tumbles(ii,2)+1),'b');
        end
        
        set(speeds_plot,'xdata',time,'ydata',speeds)
        set(timeline(3),'xdata',nan,'ydata',nan)
        
        for ii=1:size(runs,1)
            run_lines(ii + size(runs,1)) = plot(ax(3),time(runs(ii,1):runs(ii,2)),speeds(runs(ii,1):runs(ii,2)),'r');
        end
        
        for ii = 1:size(tumbles,1)
            tumble_lines(ii + size(tumbles,1)) = plot(ax(3),time(tumbles(ii,1)-1:tumbles(ii,2)+1),speeds(tumbles(ii,1)-1:tumbles(ii,2)+1),'b');
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
            
            if strcmp(ext,'avi')
                I = read(V,track(ff,2));
            elseif strcmp(ext,'tif')
                I = imread(V(1).Filename, track(ff,2));
            end
            
            set(p2,'Cdata',I)
            set(p, 'Xdata', nan, 'YData', nan)
            set(c, 'Xdata', nan, 'YData', nan)
            
            x = Geometries.Boundaries{track(ff,2)}{track(ff,3)}(:,2);
            y = Geometries.Boundaries{track(ff,2)}{track(ff,3)}(:,1);
            set(p, 'Xdata', x, 'YData', y)
            set(c, 'Xdata', Geometries.Centroids{track(ff,2)}(track(ff,3),1), 'YData', Geometries.Centroids{track(ff,2)}(track(ff,3),2))
            set(track_lines{1}, 'Xdata', track(1:ff,5)*scaling,...
                'YData', track(1:ff,6)*scaling)
            
            if ~isempty(runs) && ismember(ff,runs(:,1))
                ff0 = ff;
                ll = ll+1;
                color = 'r';
                track_lines{ll} = plot(ax(1), nan, nan, 'color',color);
            end
            if ~isempty(tumbles) && ismember(ff,tumbles(:,1))
                ff0 = ff;
                ll = ll+1;
                color = 'b';
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
            
            pause(1/FrameRate-toc(timerVal))
        end
        stop = true;
    end
end