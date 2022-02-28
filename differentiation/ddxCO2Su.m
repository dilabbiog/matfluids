%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         DIFFERENTIATION TOOLBOX                         %
%                                                                         %
% ddxCO2Su                                                                %
% Compact Difference Derivative                                           %
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
% du = ddxCO2Su(u, dx);                                                   %
% du = ddxCO2Su(u, dx, drc);                                              %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Compute the first derivative of an input array for a uniform spacing in %
% the specified direction. This function applies the classical second-    %
% order compact (or Padé) finite difference scheme described in [1]. At   %
% the boundaries, a second-order forward or backward compact difference   %
% scheme is used. The default derivative direction is 1. When a vector is %
% being differentiated, the direction is determined automatically.        %
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
% [1] Lele, S. K. (1992). Compact finite difference schemes with          %
%     spectral-like resolution. Journal of Computational Physics, 103(1), %
%     16-42.                                                              %
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
% pi] using the classical second-order compact difference scheme with a   %
% grid spacing of pi/50.                                                  %
%                                                                         %
% >> dx      = pi/50;                                                     %
% >> u       = sin(0:dx:pi).';                                            %
% >> ux_O2   = ddxCO2Su(u, dx);                                           %
% >> ux_TRUE = cos(0:dx:pi).';                                            %
% >> E2      = abs(ux_O2 - ux_TRUE);                                      %
% >> disp(max(E2(:)));                                                    %
%     8.9843e-04                                                          %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Compute the derivative of the function u = y^2*sin(x) with respect to y %
% on the interval (x,y) = ([0,pi],[0,1]) using the classical second-order %
% compact difference scheme with a grid spacing of (pi/50, 0.01).         %
%                                                                         %
% >> dx      = pi/50;                                                     %
% >> dy      = 0.01;                                                      %
% >> x       = (0:dx:pi).';                                               %
% >> y       = (0:dy:1).';                                                %
% >> [X,Y]   = ndgrid(x,y);                                               %
% >> u       = (Y.^2).*sin(X);                                            %
% >> uy_O2   = ddxCO2Su(u, dy, 2);                                        %
% >> uy_TRUE = 2*Y.*sin(X);                                               %
% >> E2      = abs(uy_O2 - uy_TRUE);                                      %
% >> disp(max(E2(:)));                                                    %
%    1.8208e-14                                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [du] = ddxCO2Su(u, dx, varargin)


%% PARSE INPUTS

% Input defaults.
default.drc = 1;

% Input checks.
check.drc = @(x) validateattributes(x,                                  ...
                 {'logical', 'numeric'},                                ...
                 {'finite', 'positive', 'scalar'});
check.dx  = @(x) validateattributes(x,                                  ...
                 {'logical', 'numeric'},                                ...
                 {'finite', 'nonzero', 'scalar'});
check.u   = @(x) validateattributes(x,                                  ...
                 {'logical', 'numeric'},                                ...
                 {'nonempty'});

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
len = sz(drc);
clear dim;

% Define the permutation order.
pm = 1:nd;
if drc > 1
    pm(drc) = [];
    pm      = [drc pm];
    u       = permute(u, pm);
end
szp = sz(pm);

% Initialize the derivative array.
du = zeros(szp);                                               %#ok<PREALL>

% Generate colon index assignments.
idx = repmat({':'},1,nd-1);

% Define the compact scheme constants for interior and boundary nodes.
inode = struct('alpha', 1/4, 'beta', 0, 'a', 3/2, 'b', 0, 'c', 0);
bnode = struct('alpha', 1, 'beta', 0, 'a', -2, 'b', 2, 'c', 0, 'd', 0);

% Define the coefficient matrix.
dl2 = [inode.beta*ones(len-3,1); bnode.beta; 0; 0];
dl1 = [inode.alpha*ones(len-2,1); bnode.alpha; 0];
du1 = [0; bnode.alpha; inode.alpha*ones(len-2,1)];
du2 = [0; 0; bnode.beta; inode.beta*ones(len-3,1)];
A   = spdiags([dl2 dl1 ones(len,1) du1 du2], -2:2, len, len);

% Initialize the constant matrix.
b = zeros(szp);

% Compute the constant matrix at the first boundary.
b(1,idx{:}) = (bnode.a*u(1,idx{:}) + bnode.b*u(2,idx{:}))/dx;

% Compute the constant matrix at the last boundary.
b(end,idx{:}) = -(bnode.a*u(end,idx{:}) + bnode.b*u(end-1,idx{:}))/dx;

% Compute the constant matrix everywhere in between.
b(2:end-1,idx{:}) = inode.a*(u(3:end,idx{:}) - u(1:end-2,idx{:}))/(2*dx);

% Compute and inverse permute the derivative.
du = A\b;
du = ipermute(du, pm);
clear drc len nd pm sz szp;


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
% Line(s) 183                                                             %
% Message(s)                                                              %
% * The preallocated value assigned to variable 'du' might be unused.     %
% Reason(s)                                                               %
% * Since we may be dealing with large n-dimensional arrays, it is best   %
%   to preallocate memory for the variable before we begin computing the  %
%   derivative.                                                           %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% CHANGE LOG                                                              %
%                                                                         %
% 2022/02/28 -- (GDL) Changed variable name: dir -> drc.                  %
% 2022/02/25 -- (GDL) Added an if statement for permuting u.              %
% 2022/02/23 -- (GDL) Separated comments for variables A and b.           %
% 2022/02/23 -- (GDL) Changed variable name n to len.                     %
% 2022/02/23 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
