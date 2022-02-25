%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         DIFFERENTIATION TOOLBOX                         %
%                                                                         %
% ddxO2Gu                                                                 %
% Finite Difference Derivative                                            %
% First derivative, second-order error, geometry, uniform spacing         %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% D�partement de g�nie m�canique                                          %
% �cole de technologie sup�rieure (�TS)                                   %
% Montr�al, Qu�bec                                                        %
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
% du = ddxO2Gu(u, dx, geom);                                              %
% du = ddxO2Gu(u, dx, geom, dir);                                         %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Compute the first derivative of an input array associated with a        %
% geometry for a uniform spacing in the specified direction. This         %
% function applies the second-order central finite difference scheme. At  %
% the boundaries, a second-order forward or backward finite difference    %
% scheme is used. In cases where there are not enough points to apply the %
% second-order scheme, this function successively switches to lower-order %
% schemes. If there are not enough points for the lower-order schemes,    %
% the derivative is simply set to zero. When a vector is being            %
% differentiated, the direction is determined automatically.              %
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
% 'geom'         LOGICAL/NUMERIC N-DIMENSIONAL ARRAY                      %
%              ~ Input array. Geometry or mask associated with the input  %
%                array. Values are 1 within the geometry and 0 outside    %
%                the geometry.                                            %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'dir'          POSITIVE INTEGER SCALAR                                  %
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
% pi] within some defined geometrical bounds using the second-order,      %
% central finite difference scheme with a grid spacing of pi/50.          %
%                                                                         %
% >> dx      = pi/50;                                                     %
% >> geom    = [ones(20,1); zeros(10,1); ones(6,1); zeros(15,1)];         %
% >> u       = geom.*(sin(0:dx:pi).');                                    %
% >> ux_O2   = ddxO2Gu(u, dx, geom);                                      %
% >> ux_TRUE = geom.*(cos(0:dx:pi).');                                    %
% >> E2      = abs(ux_O2 - ux_TRUE);                                      %
% >> disp(max(E2(:)));                                                    %
%     0.0013                                                              %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Compute the derivative of the function u = y^2*sin(x) with respect to y %
% on the interval (x,y) = ([0,pi],[0,1]) within some defined geometrical  %
% bounds using the second-order, central finite difference scheme with a  %
% grid spacing of (pi/50, 0.01).                                          %
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
% >> uy_O2             = ddxO2Gu(u, dy, geom, 2);                         %
% >> uy_TRUE           = 2*geom.*Y.*sin(X);                               %
% >> E2                = abs(uy_O2 - uy_TRUE);                            %
% >> disp(max(E2(:)));                                                    %
%    9.3259e-15                                                           %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 3                                                               %
%                                                                         %
% Repeat the previous example using sparse matrices (best practice).      %
%                                                                         %
% >> dx                = pi/50;                                           %
% >> dy                = 0.01;                                            %
% >> x                 = (0:dx:pi).';                                     %
% >> y                 = (0:dy:1).';                                      %
% >> [X,Y]             = ndgrid(x,y);                                     %
% >> geom              = zeros(size(X));                                  %
% >> geom(5:20,5:20)   = 1;                                               %
% >> geom(40:45,55:80) = 1;                                               %
% >> geom              = sparse(logical(geom));                           %
% >> u                 = geom.*(Y.^2).*sin(X);                            %
% >> uy_O2             = ddxO2Gu(u, dy, geom, 2);                         %
% >> uy_TRUE           = 2*geom.*Y.*sin(X);                               %
% >> E2                = abs(uy_O2 - uy_TRUE);                            %
% >> disp(max(E2(:)));                                                    %
%    9.3259e-15                                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [du] = ddxO2Gu(u, dx, geom, varargin)


%% PARSE INPUTS

% Input defaults.
default.dir = 1;

% Input checks.
check.dir  = @(x) validateattributes(x,                                 ...
                  {'logical', 'numeric'},                               ...
                  {'finite', 'positive', 'scalar'});
check.dx   = @(x) validateattributes(x,                                 ...
                  {'logical', 'numeric'},                               ...
                  {'finite', 'nonzero', 'scalar'});
check.geom = @(x) validateattributes(x,                                 ...
                  {'logical', 'numeric'},                               ...
                  {'nonempty', 'binary'});
check.u    = @(x) validateattributes(x,                                 ...
                  {'logical', 'numeric'},                               ...
                  {'nonempty'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'u'    , check.u                 );
addRequired ( hParser, 'dx'   , check.dx                );
addRequired ( hParser, 'geom' , check.geom              );
addOptional ( hParser, 'dir'  , default.dir , check.dir );
parse(hParser, u, dx, geom, varargin{:});
clear check default;

% Additional verifications.
nargoutchk(0,1);


%% DIFFERENTIATION

% Determine the size of the array.
sz = size(u);

% Set the derivative direction.
[dim, dir] = mathdim(u);
if dim > 1, dir = hParser.Results.dir; end

% Define the permutation order.
tmp      = 1:ndims(u);
tmp(dir) = [];
pm       = [dir tmp];
clear tmp;

% Initialize the derivative array.
du = zeros(sz(pm));

% Permute the array that is being differentiated and its geometry.
u    = permute(u, pm);
geom = permute(geom, pm);

% Find all nonzero elements in the geometry (linear index form).
idxNZ = find(geom);

% Take the difference of the modulus of the linear indices.
idxNZm              = mod(idxNZ,sz(dir));
idxNZm(idxNZm == 0) = sz(dir);
idxNZd              = [diff(idxNZm); -1];

% Determine all the indices for backward differentiation.
idx  = (idxNZd ~= 1);
idxb = idxNZ(idx);

% Determine all the indices for forward differentiation.
idxf = idxNZ([true; idx(1:end-1)]);

% Determine all the indices for centred differentiation.
idx  = ~(idx + [true; idx(1:end-1)]);
idxc = idxNZ(idx);

% Determine the number of elements in each packet of indices.
n = idxb - idxf + 1;

% Apply 2nd order forward difference at the first boundary.
c2 = (n >= 3);
du(idxf(c2)) = (-u(idxf(c2)+2) + 4*u(idxf(c2)+1) - 3*u(idxf(c2)))/(2*dx);

% Apply 1st order forward difference at the first boundary.
c1 = (n == 2);
du(idxf(c1)) = (u(idxf(c1)+1) - u(idxf(c1)))/dx;

% When n == 1, the derivative is left to be zero (the initialized value).

% Apply 2nd order backward difference at the last boundary.
du(idxb(c2)) = (u(idxb(c2)-2) - 4*u(idxb(c2)-1) + 3*u(idxb(c2)))/(2*dx);

% Apply 1st order backward difference at the last boundary.
du(idxb(c1)) = (u(idxb(c1)) - u(idxb(c1)-1))/dx;

% Apply 2nd order central difference everywhere in between.
du(idxc) = (u(idxc+1) - u(idxc-1))/(2*dx);

% Inverse permute the derivative.
du = ipermute(du, pm);


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
% 2022/02/23 -- (GDL) Removed tic-toc in example 2.                       %
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