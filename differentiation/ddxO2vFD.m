%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         DIFFERENTIATION TOOLBOX                         %
%                                                                         %
% ddxO2vFD                                                                %
% Finite Difference Scheme                                                %
% First derivative, quasi-second-order error, variable spacing            %
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
% Copyright (C) 2022 Giuseppe Di Labbio                                   %
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
% du = ddxO2vFD(u, x);                                                    %
% du = ddxO2vFD(u, x, drc);                                               %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Compute the first derivative of an input array for a variable spacing   %
% in the specified direction. This function applies a quasi-second-order  %
% central finite difference scheme. At the boundaries, a quasi-second-    %
% order forward or backward finite difference scheme is used. The default %
% derivative direction is 1. When a vector is being differentiated, the   %
% direction is determined automatically.                                  %
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
% N/A                                                                     %
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
% EXAMPLE 2                                                               %
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
% EXAMPLE 3                                                               %
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
%    5.5955e-14                                                           %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 4                                                               %
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

% Input checks.
check.u   = @(x) validateattributes(x,                                  ...
                 {'logical', 'numeric'},                                ...
                 {'nonempty'});
check.x   = @(x) validateattributes(x,                                  ...
                 {'logical', 'numeric'},                                ...
                 {'finite', 'increasing', 'vector'});
check.drc = @(x) validateattributes(x,                                  ...
                 {'logical', 'numeric'},                                ...
                 {'finite', 'positive', 'scalar'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'u'   ,               check.u   );
addRequired ( hParser, 'x'   ,               check.x   );
addOptional ( hParser, 'drc' , default.drc , check.drc );
parse(hParser, u, x, varargin{:});
clear check default;

% Additional verifications.
nargoutchk(0,1);


%% DIFFERENTIATION

% Determine the size and number of dimensions of the input array.
sz = size(u);
nd = ndims(u);

% Set the derivative direction.
[dim, drc] = mathdim(u);
if dim > 1, drc = hParser.Results.drc; end
clear dim;

% Permute the input array.
pm = 1:nd;
if drc > 1
    pm(drc) = [];
    pm      = [drc pm];
    u       = permute(u, pm);
end
szp = sz(pm);

% Ensure that x is a column vector.
if isrow(x), x = x.'; end

% Initialize the spacing array.
dx = diff(x);

% Initialize the derivative array.
du = zeros(szp);

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
