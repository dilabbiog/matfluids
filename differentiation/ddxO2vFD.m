%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         DIFFERENTIATION TOOLBOX                         %
%                                                                         %
% ddxO2vFD                                                                %
% Finite Difference Scheme                                                %
% First derivative, quasi-second-order error, arbitrary spacing           %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Département de génie mécanique                                          %
% École de technologie supérieure (ÉTS)                                   %
% Montréal, Québec                                                        %
% Canada                                                                  %
%                                                                         %
% Contributors: Giuseppe Di Labbio                                        %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Copyright (C) 2023 Giuseppe Di Labbio                                   %
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
% du = ddxO2vFD(u, dx);                                                   %
% du = ddxO2vFD(u, x);                                                    %
% du = ddxO2vFD(___, drc);                                                %
% du = ddxO2vFD(___, Name, Value);                                        %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Compute the first derivative of an input array for arbitrary spacing in %
% the specified direction. This function applies a quasi-second-order     %
% central finite difference scheme. At the boundaries, a quasi-second-    %
% order forward or backward finite difference scheme is used. The default %
% derivative direction is 1. When a vector is being differentiated, the   %
% direction is determined automatically. If the spacing is found to be    %
% uniform within some tolerance, a second-order uniform spacing scheme is %
% used.                                                                   %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2019b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
% mathdim                                                                 %
%                                                                         %
% Acknowledgments:                                                        %
% N/A                                                                     %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ======================================================================= %
% Input Arguments (Required):                                             %
% ----------------------------------------------------------------------- %
% 'u'            LOGICAL/NUMERIC N-DIMENSIONAL ARRAY                      %
%              ~ Input array. Array to be differentiated.                 %
% ----------------------------------------------------------------------- %
% 'dx'           NONZERO SCALAR                                           %
%              ~ Input scalar. Uniform spacing in the direction of the    %
%                differentiation.                                         %
% ----------------------------------------------------------------------- %
% 'x'            MONOTONICALLY INCREASING/DECREASING 1-D NUMERIC ARRAY    %
%              ~ Input vector. Vector of coordinates for the function to  %
%                be differentiated.                                       %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'drc'          POSITIVE INTEGER SCALAR                                  %
%                Default: 1 or vector direction                           %
%              ~ Input scalar. Direction along which to perform the       %
%                derivative.                                              %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% 'tol'          POSITIVE REAL SCALAR                                     %
%                Default: 1e-8                                            %
%              ~ Input scalar. Tolerance to consider spacing as uniform.  %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'du'           NUMERIC N-DIMENSIONAL ARRAY                              %
%              ~ Output array. Derivative of the input array in the       %
%                specified direction.                                     %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Compute the derivative of the function u = sin(x) on the interval [0,   %
% pi] using the second-order, central finite difference scheme with a     %
% grid spacing of pi/50.                                                  %
%                                                                         %
% >> dx      = pi/50;                                                     %
% >> u       = sin(0:dx:pi).';                                            %
% >> ux_O2   = ddxO2vFD(u, dx);                                           %
% >> ux_TRUE = cos(0:dx:pi).';                                            %
% >> E2      = abs(ux_O2 - ux_TRUE);                                      %
% >> disp(max(E2(:)));                                                    %
%     0.0013                                                              %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Compute the derivative of the function u = sin(x) on the interval [0,   %
% pi] using the quasi-second-order, central finite difference scheme with %
% a constant grid spacing of pi/50.                                       %
%                                                                         %
% >> dx      = pi/50;                                                     %
% >> x       = (0:dx:pi).';                                               %
% >> u       = sin(x);                                                    %
% >> ux_O2   = ddxO2vFD(u, x);                                            %
% >> ux_TRUE = cos(x);                                                    %
% >> E2      = abs(ux_O2 - ux_TRUE);                                      %
% >> disp(max(E2(:)));                                                    %
%     0.0013                                                              %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 3                                                               %
%                                                                         %
% Compute the derivative of the function u = y^2*sin(x) with respect to y %
% on the interval (x,y) = ([0,pi],[0,1]) using the second-order, central  %
% finite difference scheme with a grid spacing of (pi/50, 0.01).          %
%                                                                         %
% >> dx      = pi/50;                                                     %
% >> dy      = 0.01;                                                      %
% >> x       = (0:dx:pi).';                                               %
% >> y       = (0:dy:1).';                                                %
% >> [X,Y]   = ndgrid(x,y);                                               %
% >> u       = (Y.^2).*sin(X);                                            %
% >> uy_O2   = ddxO2vFD(u, dy, 2);                                        %
% >> uy_TRUE = 2*Y.*sin(X);                                               %
% >> E2      = abs(uy_O2 - uy_TRUE);                                      %
% >> disp(max(E2(:)));                                                    %
%    3.5083e-14                                                           %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 4                                                               %
%                                                                         %
% Compute the derivative of the function u = exp(-x) on the interval [2,  %
% 10] using a quasi-second-order, central finite difference scheme with a %
% grid of 51 points spaced out as a decreasing geometric series.          %
%                                                                         %
% >> x       = geospace(2, 10, 51, 'reverse');                            %
% >> u       = exp(-x);                                                   %
% >> ux_O2   = ddxO2vFD(u, x);                                            %
% >> ux_TRUE = -exp(-x);                                                  %
% >> E2      = abs(ux_O2 - ux_TRUE);                                      %
% >> disp(max(E2(:)));                                                    %
%    1.7176e-04                                                           %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 5                                                               %
%                                                                         %
% Compute the derivative of the function u = y^2*sin(x) with respect to y %
% on the interval (x,y) = ([0,pi],[0,1]) using the quasi-second-order,    %
% central finite difference scheme with a constant grid spacing of        %
% (pi/50, 0.01).                                                          %
%                                                                         %
% >> dx      = pi/50;                                                     %
% >> dy      = 0.01;                                                      %
% >> x       = (0:dx:pi).';                                               %
% >> y       = (0:dy:1).';                                                %
% >> [X,Y]   = ndgrid(x,y);                                               %
% >> u       = (Y.^2).*sin(X);                                            %
% >> uy_O2   = ddxO2vFD(u, y, 2);                                         %
% >> uy_TRUE = 2*Y.*sin(X);                                               %
% >> E2      = abs(uy_O2 - uy_TRUE);                                      %
% >> disp(max(E2(:)));                                                    %
%    5.3069e-14                                                           %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 6                                                               %
%                                                                         %
% Compute the derivative of the function u = y^2*sin(x) with respect to y %
% on the interval (x,y) = ([0,pi],[0,1]) using the quasi-second-order,    %
% central finite difference scheme with a constant grid spacing of in the %
% x direction of pi/50 and a grid of 25 points in the y direction spaced  %
% out as an increasing geometric series.                                  %
%                                                                         %
% >> dx      = pi/50;                                                     %
% >> x       = (0:dx:pi).';                                               %
% >> y       = geospace(0, 1, 25);                                        %
% >> [X,Y]   = ndgrid(x,y);                                               %
% >> u       = (Y.^2).*sin(X);                                            %
% >> uy_O2   = ddxO2vFD(u, y, 2);                                         %
% >> uy_TRUE = 2*Y.*sin(X);                                               %
% >> E2      = abs(uy_O2 - uy_TRUE);                                      %
% >> disp(max(E2(:)));                                                    %
%    6.5344e-09                                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [du] = ddxO2vFD(u, x, varargin)


%% PARSE INPUTS

% Input defaults.
default.drc = 1;
default.tol = 1e-8;

% Input checks.
if isempty(u), check.u = @(x) 1;
else,          check.u = @(x) validateattributes(x,                     ...
                              {'logical', 'numeric'},                   ...
                              {'nonempty'});
end
if isempty(x), check.x = @(x) 1;
else,          check.x = @(x) validateattributes(x,                     ...
                              {'logical', 'numeric'},                   ...
                              {'finite', 'vector'});
end
check.drc = @(x) validateattributes(x,                                  ...
                 {'logical', 'numeric'},                                ...
                 {'finite', 'positive', 'scalar'});
check.tol = @(x) validateattributes(x,                                  ...
                 {'logical', 'numeric'},                                ...
                 {'finite', 'nonnegative', 'scalar'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser , 'u'   ,               check.u   );
addRequired ( hParser , 'x'   ,               check.x   );
addOptional ( hParser , 'drc' , default.drc , check.drc );
addParameter( hParser , 'tol' , default.tol , check.tol );
parse(hParser, u, x, varargin{:});
clear check default;

% Additional verifications.
narginchk(2,5);
nargoutchk(0,1);


%% SPECIAL CASES

% If the variable u is empty/scalar, its derivative is empty/zero.
if isempty(u)
    du = [];
    return;
elseif isscalar(u)
    du = 0;
    return;
end


%% INITIALIZATIONS

% Initialize the spacing array.
if isempty(x)
    dx = 1;
elseif isscalar(x)
    dx = x;
else
    dx = diff(x);
    if isrow(dx), dx = dx.'; end
    
    % Check if coordinate is monotonically increasing/decreasing.
    sgn = sign(dx);
    if max(sgn) ~= min(sgn)
        error('coord:notMonotonic', ...
             ['Error. \nThe input coordinate must be strictly monotoni' ...
              'cally increasing or decreasing.']);
    end
    
    % Check if coordinate has uniform spacing.
    dx0 = uniquetol(dx, hParser.Results.tol);
    if isscalar(dx0), dx = dx0; end
    clear dx0;
end

% Set the derivative direction.
[dim, drc] = mathdim(u);
if dim > 1, drc = hParser.Results.drc; end
clear dim;

% Determine the size and number of dimensions of the input array.
sz = size(u);
nd = ndims(u);

% Permute the input array.
pm = 1:nd;
if drc > 1
    pm(drc) = [];
    pm      = [drc pm];
    u       = permute(u, pm);
end
szp = sz(pm);

% Initialize the derivative array.
du = zeros(szp);


%% DIFFERENTIATION (UNIFORM SPACING)

% For scalar spacing.
if isscalar(dx)
    du = ddxO2uFD(u, dx, 1);
    du = ipermute(du, pm);
    clear drc nd pm sz szp;
    return;
end


%% DIFFERENTIATION (VARIABLE SPACING)

% Generate colon index assignments.
idx = repmat({':'}, 1, nd-1);

% Apply quasi-2nd order forward difference at the first boundary.
du(1,idx{:}) = (-dx(1)^2*u(3,idx{:})                                    ...
             +  (dx(1) + dx(2))^2*u(2,idx{:})                           ...
             -   dx(2)*(2*dx(1) + dx(2))*u(1,idx{:}))                   ...
             /  (dx(1)*dx(2)*(dx(1) + dx(2)));

% Apply quasi-2nd order backward difference at the last boundary.
du(end,idx{:}) = (dx(end)^2*u(end-2,idx{:})                             ...
               - (dx(end) + dx(end-1))^2*u(end-1,idx{:})                ...
               +  dx(end-1)*(2*dx(end) + dx(end-1))*u(end,idx{:}))      ...
               / (dx(end)*dx(end-1)*(dx(end) + dx(end-1)));

% Apply quasi-2nd order central difference everywhere in between.
du(2:end-1,idx{:}) = (dx(1:end-1).^2.*u(3:end,idx{:})                   ...
                   -  dx(2:end).^2.*u(1:end-2,idx{:})                   ...
                   - (dx(1:end-1).^2 - dx(2:end).^2)                    ...
                  .* u(2:end-1,idx{:}))                                 ...
                  ./ (dx(1:end-1).*dx(2:end).*(dx(1:end-1) + dx(2:end)));

% Inverse permute the derivative.
du = ipermute(du, pm);
clear drc nd pm sz szp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  NOTES                                  %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Line(s) N/A                                                             %
% * N/A                                                                   %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% SUPPRESSED MESSAGES                                                     %
%                                                                         %
% Line(s) N/A                                                             %
% Message(s)                                                              %
% * N/A                                                                   %
% Reason(s)                                                               %
% * N/A                                                                   %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% CHANGE LOG                                                              %
%                                                                         %
% 2023/09/25 -- (GDL) Added support for uniform spacing and empty arrays. %
% 2022/03/04 -- (GDL) Modified the examples.                              %
% 2022/03/03 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
