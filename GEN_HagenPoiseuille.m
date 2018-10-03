%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           GENERATE FLUID FLOWS                          %
%                        THE HAGEN-POISEUILLE FLOW                        %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Giuseppe Di Labbio                                                      %
% Department of Mechanical, Industrial & Aerospace Engineering            %
% Concordia University Montréal, Canada                                   %
%                                                                         %
% Last Update: October 3rd, 2018 by Giuseppe Di Labbio                    %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Copyright (C) 2018 Giuseppe Di Labbio                                   %
%                                                                         %
% This program is free software: you can redistribute it and/or modify it %
% under the terms of the GNU General Public License as published by the   %
% Free Software Foundation, either version 3 of the License, or (at your  %
% option) any later version.                                              %
%                                                                         %
% This program is distributed in the hope that it will be useful, but     %
% WITHOUT ANY WARRANTY; without even the implied warranty of              %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        %
% General Public License for more details.                                %
%                                                                         %
% You should have received a copy of the GNU General Public License along %
% with this program. If not, see <https://www.gnu.org/licenses/>.         %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% SYNTAX                                                                  %
%                                                                         %
% [VEC, DRV] = GEN_HagenPoiseuille(R, dr, L, dl, t, mu, dP);              %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function generates a steady Hagen-Poiseuille flow over a grid and  %
% at the times specified by the user.                                     %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'dl'         - REAL SCALAR                                              %
%              - Discretization in the axial direction.                   %
% ----------------------------------------------------------------------- %
% 'dP'         - REAL SCALAR                                              %
%              - Pressure difference across the length of the pipe.       %
% ----------------------------------------------------------------------- %
% 'dr'         - REAL SCALAR                                              %
%              - Discretization in the radial direction.                  %
% ----------------------------------------------------------------------- %
% 'DRV'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - One-dimensional cell array of structs each containing    %
%                the derivatives of the velocity components (U and V) in  %
%                X and Y.                                                 %
% ----------------------------------------------------------------------- %
% 'L'          - REAL SCALAR                                              %
%              - Length of the pipe.                                      %
% ----------------------------------------------------------------------- %
% 'mu'         - REAL SCALAR                                              %
%              - Dynamic viscosity of the working fluid.                  %
% ----------------------------------------------------------------------- %
% 'R'          - REAL SCALAR                                              %
%              - Radius of the circular pipe.                             %
% ----------------------------------------------------------------------- %
% 't'          - REAL VECTOR                                              %
%              - Column vector of times at which to compute the double    %
%                gyre flow (t(1) = min, t(end) = max).                    %
% ----------------------------------------------------------------------- %
% 'VEC'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - One-dimensional cell array of structs each containing    %
%                information on the spatial mask (C), velocity components %
%                (U and V), and Cartesian grid (X and Y).                 %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Generate a steady Hagen-Poseuille flow for a pipe of length 10 m and    %
% radius 2 cm. Use a grid spacing of 0.02 cm in the radial direction and  %
% 0.1 m in the axial direction. Assume a dynamic viscosity of 1 cP and a  %
% pressure difference of 0.1 kPa across the length of the pipe. Make a    %
% time series from t = 0 s to t = 5 s in increments of 0.1 s.             %
%                                                                         %
% >> R  = 0.02;                                                           %
% >> dr = 0.0002;                                                         %
% >> L  = 10;                                                             %
% >> dl = 0.1;                                                            %
% >> t  = (0:0.1:5).';                                                    %
% >> mu = 0.001;                                                          %
% >> dP = 100;                                                            %
% >> [VEC, VGT] = GEN_HagenPoiseuille(R, dr, L, dl, t, mu, dP);           %
% >> quiver(VEC{1}.X(1:4:end,1:4:end), VEC{1}.Y(1:4:end,1:4:end), ...     %
%           VEC{1}.U(1:4:end,1:4:end), VEC{1}.V(1:4:end,1:4:end), ...     %
%           'ShowArrowHead', 'off');                                      %
% >> axis([0 10 -0.05 0.05]);                                             %
% >> hold on;                                                             %
% >> plot([0 10], [0.02 0.02],   'LineWidth', 3, 'Color', 'k');           %
% >> plot([0 10], [-0.02 -0.02], 'LineWidth', 3, 'Color', 'k');           %
% >> hold off;                                                            %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% N/A                                                                     %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GEN_HagenPoiseuille
function [VEC, DRV] = GEN_HagenPoiseuille(R, dr, L, dl, t, mu, dP)

% Determine the number of time steps (n).
n = length(t);

% Create the x and y vectors of the grid and determine their size.
x  = (0:dl:L).';
nx = length(x);
y  = (-R:dr:R).';
ny = length(y);

% Initialize the cell arrays to hold the velocity field and derivatives.
VEC = cell(n,1);
DRV = cell(n,1);

% Declare a struct holding data for the mask (although there is none in
% this case), the velocity components, and the Cartesian grids.
VEC{1} = struct('C', ones(ny,nx), 'U', 0, 'V', 0, 'X', 0, 'Y', 0);

% Declare a struct holding the derivatives of the velocity components.
DRV{1} = struct('UX', 0, 'UY', 0, 'VX', 0, 'VY', 0);

% Create the Cartesian grid.
[VEC{1}.X, VEC{1}.Y] = meshgrid(x,flipud(y));

% Compute the velocity field.
VEC{1}.U = (0.25/mu)*(dP/L)*(R^2 - VEC{1}.Y.^2);
VEC{1}.V = zeros(ny,nx);

% Compute the exact derivatives of the velocity field.
DRV{1}.UX = zeros(ny,nx);
DRV{1}.VX = zeros(ny,nx);
DRV{1}.UY = -(0.5/mu)*(dP/L)*VEC{1}.Y;
DRV{1}.VY = zeros(ny,nx);

% Repeat the velocity field and velocity gradient tensor values for the
% requested time series.
for k = 2:n
    VEC{k}.C = VEC{1}.C;
    VEC{k}.X = VEC{1}.X;
    VEC{k}.Y = VEC{1}.Y;
    VEC{k}.U = VEC{1}.U;
    VEC{k}.V = VEC{1}.V;
    
    DRV{k}.UX = DRV{1}.UX;
    DRV{k}.VX = DRV{1}.VX;
    DRV{k}.UY = DRV{1}.UY;
    DRV{k}.VY = DRV{1}.VY;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% SUPPRESS MESSAGES %%%%%%%%%%%%%%%%%%%%%%%%% %%

%#ok<*N/A>
% Line(s) N/A
% Message(s)
% * N/A
% Reason(s)
% * N/A

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Line(s) 134                                                             %
% * Note that the 'C' field, which should contain the mask information,   %
%   is used only to be consistent with other codes in this package.       %
%                                                                         %
% Line(s) 140                                                             %
% * Note that flipud(y) is used in constructing the grid. This is simply  %
%   in order to be consistent with other codes in this package where the  %
%   values of X increase from left to right and the values of Y increase  %
%   from bottom to top.                                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
