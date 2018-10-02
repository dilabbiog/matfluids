%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           GENERATE FLUID FLOWS                          %
%                        TIME-DEPENDENT DOUBLE GYRE                       %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Giuseppe Di Labbio                                                      %
% Department of Mechanical, Industrial & Aerospace Engineering            %
% Concordia University Montréal, Canada                                   %
%                                                                         %
% Last Update: October 2nd, 2018 by Giuseppe Di Labbio                    %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% SYNTAX                                                                  %
%                                                                         %
% [VEC, DRV] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                    %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function generates a time-dependent double gyre flow over a grid   %
% and at the times specified by the user.                                 %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'A'          - REAL SCALAR                                              %
%              - Amplitude or velocity scale.                             %
% ----------------------------------------------------------------------- %
% 'DRV'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - One-dimensional cell array of structs each containing    %
%                the derivatives of the velocity components (U and V) in  %
%                X and Y.                                                 %
% ----------------------------------------------------------------------- %
% 'epsn'       - REAL SCALAR                                              %
%              - An approximation of the motion amplitude of the gyres.   %
% ----------------------------------------------------------------------- %
% 'omga'       - REAL SCALAR                                              %
%              - Radial frequency of oscillation (period = 2*pi/omga).    %
% ----------------------------------------------------------------------- %
% 't'          - REAL VECTOR                                              %
%              - Column vector of times at which to compute the double    %
%                gyre flow (t(1) = min, t(end) = max).                    %
% ----------------------------------------------------------------------- %
% 'x'          - REAL VECTOR                                              %
%              - Column vector of x coordinates over which to compute the %
%                double gyre flow (x(1) = min, x(end) = max).             %
% ----------------------------------------------------------------------- %
% 'y'          - REAL VECTOR                                              %
%              - Column vector of y coordinates over which to compute the %
%                double gyre flow (y(1) = min, y(end) = max).             %
% ----------------------------------------------------------------------- %
% 'VEC'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - One-dimensional cell array of structs each containing    %
%                information on the spatial mask (C), velocity components %
%                (U and V), and Cartesian grid (X and Y).                 %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Generate a time-dependent double gyre on the domain (x,y) = [0,2]x[0,1] %
% with a constant grid spacing of 0.01 over the time interval [0,20] with %
% time-step size 0.1. Use A = 0.1, epsilon = 0.25, and omega = 2*pi/10.   %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A    = 0.1;                                                          %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> quiver(VEC{1}.X(1:4:end,1:4:end), VEC{1}.Y(1:4:end,1:4:end), ...     %
%           VEC{1}.U(1:4:end,1:4:end), VEC{1}.V(1:4:end,1:4:end));        %
% >> axis([0 2 0 1]);                                                     %
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

%% GEN_DoubleGyre
function [VEC, DRV] = GEN_DoubleGyre(x, y, t, A, epsn, omga)

% Determine the number of time steps (n) and the number of data points in
% the x (nx) and y (ny) directions.
n  = length(t);
nx = length(x);
ny = length(y);

% Initialize the cell arrays to hold the velocity field and derivatives.
VEC = cell(n,1);
DRV = cell(n,1);
for k = 1:n
    
    % Declare a struct holding data for the mask (although there is none in
    % this case), the velocity components, and the Cartesian grids.
    VEC{k} = struct('C', ones(ny,nx), 'U', 0, 'V', 0, 'X', 0, 'Y', 0);
    
    % Declare a struct holding the derivatives of the velocity components.
    DRV{k} = struct('UX', 0, 'UY', 0, 'VX', 0, 'VY', 0);
    
    % Create the Cartesian grid.
    [VEC{k}.X,VEC{k}.Y] = meshgrid(x,flipud(y));
    
    % Compute the velocity field.
    a = epsn*sin(omga*t(k));
    b = 1 - 2*epsn*sin(omga*t(k));
    f = a*VEC{k}.X.^2 + b*VEC{k}.X;
    VEC{k}.U =  -A*pi*sin(pi*f).*cos(pi*VEC{k}.Y);
    VEC{k}.V = A*pi*(2*a*VEC{k}.X + b).*cos(pi*f).*sin(pi*VEC{k}.Y);
    
    % Compute the exact derivatives of the velocity field.
    DRV{k}.UX = -A*pi^2*(2*a*VEC{k}.X + b).*cos(pi*f).*cos(pi*VEC{k}.Y);
    DRV{k}.VX = A*pi*(2*a*cos(pi*f) - pi*(2*a*VEC{k}.X + b).^2 ...
                .*sin(pi*f)).*sin(pi*VEC{k}.Y);
    DRV{k}.UY = A*pi^2*sin(pi*f).*sin(pi*VEC{k}.Y);
    DRV{k}.VY = A*pi^2*(2*a*VEC{k}.X + b).*cos(pi*f).*cos(pi*VEC{k}.Y);
    
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
% Line(s) 106                                                             %
% * Note that the 'C' field, which should contain the mask information,   %
%   is used only to be consistent with other codes in this package.       %
%                                                                         %
% Line(s) 112                                                             %
% * Note that flipud(y) is used in constructing the grid. This is simply  %
%   in order to be consistent with other codes in this package where the  %
%   values of X increase from left to right and the values of Y increase  %
%   from bottom to top.                                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
