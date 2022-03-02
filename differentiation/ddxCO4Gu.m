%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         DIFFERENTIATION TOOLBOX                         %
%                                                                         %
% ddxCO4Gu                                                                %
% Compact Difference Derivative                                           %
% First derivative, fourth-order error, geometry, uniform spacing         %
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
% du = ddxCO4Gu(u, dx, geom);                                             %
% du = ddxCO4Gu(u, dx, geom, drc);                                        %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Compute the first derivative of an input array associated with a        %
% geometry for a uniform spacing in the specified direction. This         %
% function applies the classical fourth-order compact (or Padé) finite    %
% difference scheme described in [1]. At the boundaries, a fourth-order   %
% forward or backward compact difference scheme is used. In cases where   %
% there are not enough points to apply the fourth-order scheme, this      %
% function successively switches to lower-order schemes. If there are not %
% enough points for the lower-order schemes, the derivative is simply set %
% to zero. The default derivative direction is 1. When a vector is being  %
% differentiated, the direction is determined automatically.              %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2019b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
% idNodeidx                                                               %
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
% ----------------------------------------------------------------------- %
% 'geom'         LOGICAL/NUMERIC N-DIMENSIONAL ARRAY                      %
%              ~ Input array. Geometry or mask associated with the input  %
%                array. Values are 1 within the geometry and 0 outside    %
%                the geometry.                                            %
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
% pi] within some defined geometrical bounds using the classical fourth-  %
% order compact difference scheme with a grid spacing of pi/50.           %
%                                                                         %
% >> dx      = pi/50;                                                     %
% >> geom    = [ones(20,1); zeros(10,1); ones(6,1); zeros(15,1)];         %
% >> u       = geom.*(sin(0:dx:pi).');                                    %
% >> ux_O4   = ddxCO4Gu(u, dx, geom);                                     %
% >> ux_TRUE = geom.*(cos(0:dx:pi).');                                    %
% >> E4      = abs(ux_O4 - ux_TRUE);                                      %
% >> disp(max(E4(:)));                                                    %
%    2.2836e-06                                                           %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Compute the derivative of the function u = y^2*sin(x) with respect to y %
% on the interval (x,y) = ([0,pi],[0,1]) within some defined geometrical  %
% bounds using the classical fourth-order compact difference scheme with  %
% a grid spacing of (pi/50, 0.01).                                        %
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
% >> uy_O4             = ddxCO4Gu(u, dy, geom, 2);                        %
% >> uy_TRUE           = 2*geom.*Y.*sin(X);                               %
% >> E4                = abs(uy_O4 - uy_TRUE);                            %
% >> disp(max(E4(:)));                                                    %
%    4.9738e-14                                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [du] = ddxCO4Gu(u, dx, geom, varargin)


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
du = zeros(szp);

% Determine the leading, interior and trailing node indices.
[idxf, idxc, idxb, ~] = idNodeidx(geom);

% Determine the number of elements in each packet of indices.
n    = idxb - idxf + 1;
npck = length(n);
ncum = cumsum(n);
nsum = ncum(end);

% Define the compact scheme constants for interior and boundary nodes.
inode4 = struct('alpha', 1/4, 'beta', 0, 'a', 3/2, 'b', 0, 'c', 0);
bnode4 = struct('alpha', 3, 'beta', 0,                                  ...
                'a', -17/6, 'b', 3/2, 'c', 3/2, 'd', -1/6);
bnode3 = struct('alpha', 2, 'beta', 0,                                  ...
                'a', -5/2, 'b', 2, 'c', 1/2, 'd', 0);
bnode2 = struct('alpha', 1, 'beta', 0, 'a', -2, 'b', 2, 'c', 0, 'd', 0);

% Define the coefficient matrix.
dl2 = zeros(nsum,1);
dl1 = zeros(nsum,1);
du1 = zeros(nsum,1);
du2 = zeros(nsum,1);
idx = 1;
for k = 1:npck
    dl2(idx:idx+n(k)-3)   = [inode4.beta*ones(n(k)-3,1); bnode4.beta];
    dl1(idx:idx+n(k)-2)   = [inode4.alpha*ones(n(k)-2,1); bnode4.alpha];
    du1(idx+1:idx+n(k)-1) = [bnode4.alpha; inode4.alpha*ones(n(k)-2,1)];
    du2(idx+2:idx+n(k)-1) = [bnode4.beta; inode4.beta*ones(n(k)-3,1)];
    idx                   = idx + n(k);
end
A = spdiags([dl2 dl1 ones(nsum,1) du1 du2], -2:2, nsum, nsum);
clear k idx dl2 dl1 du1 du2;

% Initialize the constant matrix.
b = zeros(nsum,1);

% Classify the indices for the constant matrix.
bidxb = ncum;
bidxf = [1; bidxb(1:end-1)+1];
bidxc = cell(npck,1);
for k = 1:npck, bidxc{k} = (bidxf(k)+1:bidxb(k)-1).'; end
bidxc = cell2mat(bidxc);
clear k;

% Establish conditions for boundary nodes.
bC4 = (n >  4);
bC3 = (n == 4);
bC2 = (n == 3);
bC1 = (n == 2);

% bC4: 4th-order compact scheme at these boundary nodes.
b(bidxf(bC4)) =  (bnode4.a*u(idxf(bC4))   + bnode4.b*u(idxf(bC4)+1)     ...
              +   bnode4.c*u(idxf(bC4)+2) + bnode4.d*u(idxf(bC4)+3))/dx;
b(bidxb(bC4)) = -(bnode4.a*u(idxb(bC4))   + bnode4.b*u(idxb(bC4)-1)     ...
              +   bnode4.c*u(idxb(bC4)-2) + bnode4.d*u(idxb(bC4)-3))/dx;
clear bC4 bnode4;

% bC3: 3rd-order compact scheme at these boundary nodes.
b(bidxf(bC3)) =  (bnode3.a*u(idxf(bC3)) + bnode3.b*u(idxf(bC3)+1)       ...
              +   bnode3.c*u(idxf(bC3)+2))/dx;
idx           = sub2ind(nsum*[1 1], bidxf(bC3), bidxf(bC3)+1);
A(idx)        = bnode3.alpha;
b(bidxb(bC3)) = -(bnode3.a*u(idxb(bC3)) + bnode3.b*u(idxb(bC3)-1)       ...
              +   bnode3.c*u(idxb(bC3)-2))/dx;
idx           = sub2ind(nsum*[1 1], bidxb(bC3), bidxb(bC3)-1);
A(idx)        = bnode3.alpha;
clear bC3 bnode3;

% bC2: 2nd-order compact scheme at these boundary nodes.
b(bidxf(bC2)) =  (bnode2.a*u(idxf(bC2)) + bnode2.b*u(idxf(bC2)+1))/dx;
idx           = sub2ind(nsum*[1 1], bidxf(bC2), bidxf(bC2)+1);
A(idx)        = bnode2.alpha;
b(bidxb(bC2)) = -(bnode2.a*u(idxb(bC2)) + bnode2.b*u(idxb(bC2)-1))/dx;
idx           = sub2ind(nsum*[1 1], bidxb(bC2), bidxb(bC2)-1);
A(idx)        = bnode2.alpha;
clear bC2 bnode2;

% bC1: 1st-order explicit scheme at these boundary nodes.
b(bidxf(bC1)) = (u(idxf(bC1)+1) - u(idxf(bC1)))/dx;
idx           = sub2ind(nsum*[1 1], bidxf(bC1), bidxf(bC1)+1);
A(idx)        = 0;
b(bidxb(bC1)) = (u(idxb(bC1)) - u(idxb(bC1)-1))/dx;
idx           = sub2ind(nsum*[1 1], bidxb(bC1), bidxb(bC1)-1);
A(idx)        = 0;
clear bC1;

% iC4: 4th-order compact scheme at these interior nodes.
b(bidxc) = inode4.a*(u(idxc+1) - u(idxc-1))/(2*dx);
clear inode4 idx n ncum npck nsum;

% Compute the derivative in sparse form.
spdu = A\b;
clear A b;

% Assign the indices to the derivative in full form.
du(idxb) = spdu(bidxb);
du(idxc) = spdu(bidxc);
du(idxf) = spdu(bidxf);
clear bidxb bidxc bidxf idxb idxc idxf spdu;

% Inverse permute the derivative.
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
% 2022/03/02 -- (GDL) Adjusted lower-order boundary schemes.              %
% 2022/03/02 -- (GDL) Changed function name (it is fourth order).         %
% 2022/03/02 -- (GDL) Added missing support for singular nodes.           %
% 2022/02/28 -- (GDL) Added lower-order scheme for degenerate case.       %
% 2022/02/28 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
