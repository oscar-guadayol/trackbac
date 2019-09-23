function fig = runs_and_tumbles_figure(Run_speeds, Run_lengths, Run_complete,...
    Tumble_times, Tumble_angles, FrameRate, win_length, figure_export, figname)
         
% Creates a figure with the frequency  distribution of runs and tumbles in
% the input tracks.
%
% Input variables:
%
%       Run_speeds is an MX5 matrix, where M is the number of complete runs.
%           The first column is the track ID, the 2nd column is the mean
%           speed, 3rd is the standard deviation, 4th is the minimum speed
%           and 5th is the maximum speed.
%
%       Run_lengths is an MX2 matrix, where M is the number of complete runs.
%           The first column is the track ID, and the second column is the
%           length of the run in seconds.
%
%       Tumble_times is vector with the length (in seconds) of all
%           detected tumbles.
%
%       Tumble_angles is vector with the angles (in radians) between
%           successive runs.
%
%       FrameRate is the sampling rate of the video in frames per second.
%
%       win_length is the window length (in frames) of the smoothing
%           algorithm.
%
%       figure_export is a logical value. If true it will save the figure
%           as eps file.
%
%       figname is the full path name of the exported figure.
%
% Output variables:
%
%       fig is the handle for the output figure.%
%
% Requires: 
%
% Called by: Runs_gui, Runs
%
% Ex.:
%    fig = runs_and_tumbles_figure(Run_speeds, Run_lengths, Tumble_times,...
%           Tumble_angles, FrameRate, win_length, true, figname);
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
%  phase-contrast microscope.

fig = figure('position', [1 27 1920 950]);
ax(1) = axes('position',[0.1 0.6 0.37 0.37]);
ax(2) = axes('position',[0.55 0.6 0.37 0.37]);
ax(3) = axes('position',[0.1 0.1 0.37 0.37]);
ax(4) = axes('position',[0.55 0.1 0.37 0.37]);

bin{1} = 0:win_length/FrameRate:2;
bin{2} = 0:1:ceil(max(Run_speeds(:,2)));
bin{3} = 0:win_length/FrameRate:2;

h{1} = histc(Run_lengths(logical(Run_complete(:,2)),2), bin{1});
h{2} = histc(Run_speeds(:,2),bin{2});
% h{2} = histc(Run_speeds(logical(Run_complete(:,2)),2),bin{2});
h{3} = histc(Tumble_times, bin{3});

bar(ax(1), bin{1}+0.1,h{1}./(sum(h{1}))*100,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0])
bar(ax(2), bin{2}+0.1,h{2}./(sum(h{2}))*100,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0])
bar(ax(3), bin{3}+0.1,h{3}./(sum(h{3}))*100,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0])

[tout,rout] = rose(Tumble_angles,24);
polar(ax(4),0,30,'-k') % sets the axis limit
hold(ax(4),'on')
h{4} = polar(ax(4),tout, rout./(sum(rout)/2)*100,'k');

set(fig,'position', [250   148   600   475])

%% fits exponential distributions
P{1} = fitdist(Run_lengths(logical(Run_complete(:,2)),2),'exponential');
hold(ax(1),'on')
plot(ax(1),bin{1}(1:end-1)+0.1,diff(expcdf(bin{1}, P{1}.mu))*100,'k')

P{2} = fitdist(Run_speeds(:,2),'lognormal');
hold(ax(2),'on')
plot(ax(2),bin{2}(1:end-1)+0.5,diff(cdf('lognormal',bin{2},P{2}.mu,P{2}.sigma))*100,'k')

P{3} = fitdist(Tumble_times(~isnan(Tumble_times)),'exponential');
hold(ax(3),'on')
plot(ax(3),bin{3}(1:end-1)+0.1,diff(expcdf(bin{3}, P{3}.mu))*100,'k')

xlabel(ax(1), 'Run time (s)', 'fontsize', 10);
ylabel(ax(1), 'frequency (%)', 'fontsize', 10);
ylabel(ax(1), 'frequency (%)', 'fontsize', 10);
xlabel(ax(2), 'Run speed (\mum s^{-1})', 'fontsize', 10);
xlabel(ax(3), 'Tumble time (s)', 'fontsize', 10);
ylabel(ax(3), 'frequency (%)', 'fontsize', 10);
ylabel(ax(4), 'frequency (%)', 'fontsize', 10);
ylabel(ax(4), 'Tumble angle (^o)', 'fontsize', 10);

ax(1).XLim = [min(Run_lengths(logical(Run_complete(:,2)),2)) 4];
ax(2).XLim(1) = 0;
view(ax(4),90, -90)

text(0.5, 0.5,...
    ['N = ', num2str(sum(~isnan(Run_lengths(logical(Run_complete(:,2)),2)))), newline,...
    '\mu = ', num2str(nanmean(Run_lengths(logical(Run_complete(:,2)),2)),2), 's', newline,...
    '\sigma= ', num2str(nanstd(Run_lengths(logical(Run_complete(:,2)),2)),2), 's'],...
    'Units', 'normalized',...
    'fontsize',8,...
    'parent',ax(1))
text(0.5, 0.5,...
    ['N = ', num2str(sum(~isnan(Run_speeds(:,2)))), newline,...
    '\mu = ', num2str(nanmean(Run_speeds(:,2)),4), '\mum s^{-1}', newline,...
    '\sigma= ', num2str(nanstd(Run_speeds(:,2)),3), '\mum s^{-1}'],...
    'Units', 'normalized',...
    'fontsize',8,...
    'parent',ax(2))  
text(0.5, 0.5,...
    ['N = ', num2str(size(Tumble_angles(~isnan(Tumble_angles)),1)), newline,...
     '\mu = ', num2str(nanmean(Tumble_times),2), 's', newline,...
    '\sigma= ', num2str(nanstd(Tumble_times),2), 's'],...
    'Units', 'normalized',...
    'fontsize',8,...
    'parent',ax(3))
text(0.88, 0.90,'A',...
    'Units', 'normalized',...
    'fontsize',20,...
    'parent',ax(1))
text(0.88, 0.90,'B',...
    'Units', 'normalized',...
    'fontsize',20,...
    'parent',ax(2))
text(0.88, 0.90,'C',...
    'Units', 'normalized',...
    'fontsize',20,...
    'parent',ax(3))   
text(0.88, 0.90,'D',...
    'Units', 'normalized',...
    'fontsize',20,...
    'parent',ax(4))   

%% figure export
if figure_export
    load plot_style
    plot_style.Color= 'rgb';
    plot_style.Width = 21.40;
    plot_style.Height = 16.84;
    plot_style.Resolution = 600;
    plot_style.FontName = 'Times';
    plot_style.FixedFontSize = '12';
    plot_style.Format = 'eps';
     
    if ~exist('figname','var')==1 || isempty(figname)
        [fn,fp] = uiputfile('*.eps','Save figure');
    end
    figname = [fp, fn];
    [fp,fn] = fileparts(figname);

    hgexport(fig,[fp,fn,'.eps'],plot_style)
    
    plot_style.Format = 'jpeg';
    hgexport(fig,[fp,fn,'.jpg'],plot_style)
    
end