%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         DIFFERENTIATION TOOLBOX                         %
%                                                                         %
% ddxO4uCRg                                                               %
% Hybrid Compact-Richardson Scheme                                        %
% First derivative, fourth-order error, uniform spacing, geometry         %
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
% du = ddxO4uCRg(u, dx, geom);                                            %
% du = ddxO4uCRg(u, dx, geom, drc);                                       %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Compute the first derivative of an input array associated with a        %
% geometry for a uniform spacing in the specified direction. This         %
% function applies the noise-optimized fourth-order hybrid compact-       %
% Richardson scheme [1]. The scheme makes use of the classical Padé       %
% (fourth-order compact) finite difference scheme [2] and Richardson      %
% extrapolation [3]. At the boundaries, a fourth-order forward or         %
% backward compact difference scheme is used [2]. In cases where there    %
% are not enough points to apply the fourth-order scheme, this function   %
% successively switches to lower-order schemes. If there are not enough   %
% points for the lower-order schemes, the derivative is simply set to     %
% zero. The default derivative direction is 1. When a vector is being     %
% differentiated, the direction is determined automatically.              %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2019b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
% ddxO4uCDg                                                               %
% mathdim                                                                 %
%                                                                         %
% Acknowledgments:                                                        %
% N/A                                                                     %
%                                                                         %
% References:                                                             %
% [1] Etebari, A., & Vlachos, P. P. (2005). Improvements on the accuracy  %
%     of derivative estimation from DPIV velocity measurements.           %
%     Experiments in Fluids, 39(6), 1040-1050.                            %
% [2] Lele, S. K. (1992). Compact finite difference schemes with          %
%     spectral-like resolution. Journal of Computational Physics, 103(1), %
%     16-42.                                                              %
% [3] Foucaut, J. M., & Stanislas, M. (2002). Some considerations on the  %
%     accuracy and frequency response of some derivative filters applied  %
%     to particle image velocimetry vector fields. Measurement Science    %
%     and Technology, 13(7), 1058-1071.                                   %
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
% pi] within some defined geometrical bounds using the noise-optimized    %
% fourth-order hybrid compact-Richardson scheme with a grid spacing of    %
% pi/50.                                                                  %
%                                                                         %
% >> dx      = pi/50;                                                     %
% >> geom    = [ones(20,1); zeros(10,1); ones(6,1); zeros(15,1)];         %
% >> u       = geom.*(sin(0:dx:pi).');                                    %
% >> ux_O4   = ddxO4uCRg(u, dx, geom);                                    %
% >> ux_TRUE = geom.*(cos(0:dx:pi).');                                    %
% >> E4      = abs(ux_O4 - ux_TRUE);                                      %
% >> disp(max(E4(:)));                                                    %
%     0.0258                                                              %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Compute the derivative of the function u = y^2*sin(x) with respect to y %
% on the interval (x,y) = ([0,pi],[0,1]) within some defined geometrical  %
% bounds using the noise-optimized fourth-order hybrid compact-Richardson %
% scheme with a grid spacing of (pi/50, 0.01).                            %
%                                                                         %
% >> dx                = pi/50;                                           %
% >> dy                = 0.01;                                            %
% >> x                 = (0:dx:pi).';                                     %
% >> y                 = (0:dy:1).';                                      %
% >> [X,Y]             = ndgrid(x,y);                                     %
% >> geom              = zeros(size(X));                                  %
% >> geom(5:20,5:20)   = 1;                                               %
% >> geom(40:45,55:80) = 1;                                               %
% >> u                 = geom.*(Y.^2).*sin(X);                            %
% >> uy_O4             = ddxO4uCRg(u, dy, geom, 2);                       %
% >> uy_TRUE           = 2*geom.*Y.*sin(X);                               %
% >> E4                = abs(uy_O4 - uy_TRUE);                            %
% >> disp(max(E4(:)));                                                    %
%    5.1403e-14                                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [du] = ddxO4uCRg(u, dx, geom, varargin)


%% PARSE INPUTS

% Input defaults.
default.drc = 1;

% Input checks.
check.u    = @(x) validateattributes(x,                                 ...
                  {'logical', 'numeric'},                               ...
                  {'nonempty'});
check.dx   = @(x) validateattributes(x,                                 ...
                  {'logical', 'numeric'},                               ...
                  {'finite', 'nonzero', 'scalar'});
check.geom = @(x) validateattributes(x,                                 ...
                  {'logical', 'numeric'},                               ...
                  {'nonempty', 'binary'});
check.drc  = @(x) validateattributes(x,                                 ...
                  {'logical', 'numeric'},                               ...
                  {'finite', 'positive', 'scalar'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'u'    , check.u                 );
addRequired ( hParser, 'dx'   , check.dx                );
addRequired ( hParser, 'geom' , check.geom              );
addOptional ( hParser, 'drc'  , default.drc , check.drc );
parse(hParser, u, dx, geom, varargin{:});
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

% Define the permutation order.
pm = 1:nd;
if drc > 1
    pm(drc) = [];
    pm      = [drc pm];
    u       = permute(u, pm);
    geom    = permute(geom, pm);
end
szp = sz(pm);

% Initialize the derivative array.
du = zeros(szp);                                               %#ok<PREALL>

% Generate colon index assignments.
idx = repmat({':'}, 1, nd-1);

% Define the Richardson extrapolation constants.
A  = 1239;
A1 = 272;
A2 = 1036;
A4 = -69;

% Compute the derivative.
du                  = A1*ddxO4uCDg(u, dx, geom);
du(1:2:end, idx{:}) = du(1:2:end, idx{:})                               ...
                    + A2*ddxO4uCDg(u(1:2:end, idx{:}), 2*dx,            ...
                                geom(1:2:end, idx{:}));
du(2:2:end, idx{:}) = du(2:2:end, idx{:})                               ...
                    + A2*ddxO4uCDg(u(2:2:end, idx{:}), 2*dx,            ...
                                geom(2:2:end, idx{:}));
du(1:4:end, idx{:}) = du(1:4:end, idx{:})                               ...
                    + A4*ddxO4uCDg(u(1:4:end, idx{:}), 4*dx,            ...
                                geom(1:4:end, idx{:}));
du(2:4:end, idx{:}) = du(2:4:end, idx{:})                               ...
                    + A4*ddxO4uCDg(u(2:4:end, idx{:}), 4*dx,            ...
                                geom(2:4:end, idx{:}));
du(3:4:end, idx{:}) = du(3:4:end, idx{:})                               ...
                    + A4*ddxO4uCDg(u(3:4:end, idx{:}), 4*dx,            ...
                                geom(3:4:end, idx{:}));
du(4:4:end, idx{:}) = du(4:4:end, idx{:})                               ...
                    + A4*ddxO4uCDg(u(4:4:end, idx{:}), 4*dx,            ...
                                geom(4:4:end, idx{:}));
du                  = (1/A)*du;
clear idx A A1 A2 A4;

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
% Line(s) 207                                                             %
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
% 2022/03/03 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
