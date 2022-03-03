%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         DIFFERENTIATION TOOLBOX                         %
%                                                                         %
% ddxO2u                                                                  %
% Finite Difference Derivative                                            %
% First derivative, second-order error, uniform spacing                   %
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
% du = ddxO2u(u, dx);                                                     %
% du = ddxO2u(u, dx, drc);                                                %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Compute the first derivative of an input array for a uniform spacing in %
% the specified direction. This function applies the second-order central %
% finite difference scheme. At the boundaries, a second-order forward or  %
% backward finite difference scheme is used. The default derivative       %
% direction is 1. When a vector is being differentiated, the direction is %
% determined automatically.                                               %
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
% pi] using the second-order, central finite difference scheme with a     %
% grid spacing of pi/50.                                                  %
%                                                                         %
% >> dx      = pi/50;                                                     %
% >> u       = sin(0:dx:pi).';                                            %
% >> ux_O2   = ddxO2u(u, dx);                                             %
% >> ux_TRUE = cos(0:dx:pi).';                                            %
% >> E2      = abs(ux_O2 - ux_TRUE);                                      %
% >> disp(max(E2(:)));                                                    %
%     0.0013                                                              %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
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
% >> uy_O2   = ddxO2u(u, dy, 2);                                          %
% >> uy_TRUE = 2*Y.*sin(X);                                               %
% >> E2      = abs(uy_O2 - uy_TRUE);                                      %
% >> disp(max(E2(:)));                                                    %
%    3.5083e-14                                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [du] = ddxO2u(u, dx, varargin)


%% PARSE INPUTS

% Input defaults.
default.drc = 1;

% Input checks.
check.u   = @(x) validateattributes(x,                                  ...
                 {'logical', 'numeric'},                                ...
                 {'nonempty'});
check.dx  = @(x) validateattributes(x,                                  ...
                 {'logical', 'numeric'},                                ...
                 {'finite', 'nonzero', 'scalar'});
check.drc = @(x) validateattributes(x,                                  ...
                 {'logical', 'numeric'},                                ...
                 {'finite', 'positive', 'scalar'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'u'   , check.u                 );
addRequired ( hParser, 'dx'  , check.dx                );
addOptional ( hParser, 'drc' , default.drc , check.drc );
parse(hParser, u, dx, varargin{:});
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

% Initialize the derivative array.
du = zeros(szp);

% Generate colon index assignments.
idx = repmat({':'}, 1, nd-1);

% Apply 2nd order forward difference at the first boundary.
du(1,idx{:}) = (-u(3,idx{:}) + 4*u(2,idx{:}) - 3*u(1,idx{:}))/(2*dx);

% Apply 2nd order backward difference at the last boundary.
du(end,idx{:}) = (u(end-2,idx{:}) - 4*u(end-1,idx{:})                   ...
               + 3*u(end,idx{:}))/(2*dx);

% Apply 2nd order central difference everywhere in between.
du(2:end-1,idx{:}) = (u(3:end,idx{:}) - u(1:end-2,idx{:}))/(2*dx);

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
% 2022/03/03 -- (GDL) Changed function name (removed 'S' identifier).     %
% 2022/02/28 -- (GDL) Changed input check order.                          %
% 2022/02/28 -- (GDL) Changed variable name: dir -> drc.                  %
% 2022/02/25 -- (GDL) Added an if statement for permuting u.              %
% 2022/02/23 -- (GDL) Removed message suppression in file, prefer line.   %
% 2022/02/22 -- (GDL) Formatted the code.                                 %
% 2021/07/30 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
