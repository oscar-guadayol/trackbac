function [ft,fr, Dt, Dr]  = theoretical_friction_coefficients(a, b, c, T, S, eta)
% Calculates translational and rotational friction coefficients, as well as
% translational and rotational diffusivities  using equations in appendix 6
% in Dusenbery 2009
% 
% input variables:
%
%       a, b, c are the three semiaxes of the ellipsoid (in m)
%       T is temperature in Celsius. Default = 33
%       S is salinity
%       eta is the dynamic viscosity in kg/m/s. Default =  7.4889e-04
%
% output variables:
%
%       ft is a (4X1) vector of the translation friction coefficients
%           parallel to each of the axes and the effective frictional
%           coefficient calculated as the harmonic mean of the
%           coefficients of each of the three orthogonal angles (Dusenbery
%           2009, Appendix 6).
%       fr is a (3X1) vector of the rotational friction coefficients
%           about each of the axes
%       Dt is a (4X1) vector of the translational diffusivities plus the
%           diffusion averaged over all orientations (Dusenbery
%           2009, Appendix 6).
%       Dr is a (3X1) vector of the translational diffusivities
%
%
% Requires: SW_viscosity from Thermophysical properties of seawater toolbox
%        (http://web.mit.edu/seawater/)
%
% Ex.:  [ft, fr, Dt, Dr]  = theoretical_friction_coefficients(a,b,b)
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

%% constants
if nargin<4 || isempty(T)
    T = 33; % temperature
end
if nargin<5 || isempty(S)
    S = 0;
%    eta = SW_Viscosity(T,'C',0,'ppt');
end
if nargin<6 || isempty(eta)
    eta = SW_Viscosity(T,'C',S,'ppt'); % kg/m-s
end

a = a(:);
b = b(:);
c = c(:);

k_b = 1.38e-23; % Boltzmann constant (J/K) or (m^2 kg s^-2 K^-1)
K = T+273.15; % temperature in Kelvins

%% Integrals for flow around ellipsoids (Dusenberry 2009, appendix 6)

Sint = @(x) ((a.^2 + x).*(b.^2 + x).*(c.^2 + x)).^(-0.5);
S = integral(Sint,0,Inf,'ArrayValued',true);

r = [a, b, c];
Gint = @(x) ((repmat(a,1,3).^2 + x).*(repmat(b,1,3).^2 + x).*(repmat(c,1,3).^2 + x)).^(-0.5) .* 1./(r.^2 + x);
G_E = r.^2 .* integral(Gint,0,Inf,'ArrayValued',true);

H_Ea = (G_E(:,2)+G_E(:,3))./(b.^2+c.^2); % Dusenberry 1998
H_Eb = (G_E(:,1)+G_E(:,3))./(a.^2+c.^2); % Dusenberry 1998
H_Ec = (G_E(:,1)+G_E(:,2))./(a.^2+b.^2); % Dusenberry 1998

%% transational
ft(:,1) = 16*pi*eta./(S+G_E(:,1)); % parallel to axis a eq A6.10
% ft(1) = 16*pi*eta*(a^2-b^2)/((2*a^2-b^2)*S-2*a); % eq A6.11

ft(:,2) = 16*pi*eta./(S+G_E(:,2)); %  parallel to axis b, eq A6.10
% ft(2) = 32*pi*eta*(a^2-b^2)/((2*a^2-3*b^2)*S+2*a); % eq A6.12

ft(:,3) = 16*pi*eta./(S+G_E(:,3)); %  parallel to axis c, eq A6.10
% ft(3) = 32*pi*eta*(a^2-b^2)/((2*a^2-3*b^2)*S+2*a); % eq A6.12

ft(:,4) = 3./sum(1./ft(:,1:3),2); % harmonic mean (Dusenbery 2009 appendix 6)

% from fig 4.5 in Berg 1983
% ft(1) = 4*pi*eta*a/(log(2*a/b)-1/2);
% ft(2) = 8*pi*eta*a/(log(2*a/b)+1/2);

%% rotational
fr(:,1)= 16*pi*eta/3./H_Ea; % Friction constant for rotation about the major semiaxis
fr(:,2)= 16*pi*eta/3./H_Eb; % Friction constant for rotation about the minor semiaxes
fr(:,3)= 16*pi*eta/3./H_Ec; % Friction constant for rotation about the minor semiaxes

% from fig 6.7 in Berg 1983
% fr(1) = 16/3*pi*eta*a*b^2; %kg m^2 sec^-1
% fr(2) = 8*pi*eta*a^3/3/(log(2*a/b)-1/2); %kg m^2 sec^-1

%% Difusion coefficients
Dt =  k_b*K./ft;
Dr =  k_b*K./fr;