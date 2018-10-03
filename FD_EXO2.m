%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            FINITE DIFFERENCES                           %
%                        EXPLICIT 2ND ORDER SCHEME                        %
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
% ux = FD_EXO2(u, dx);                                                    %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function applies the second-order explicit finite difference       %
% scheme to compute the derivative of a one-dimensional array. The        %
% boundary node scheme is taken as second-order as well. In cases where   %
% there are not enough points to apply the total second-order scheme,     %
% this function switches to a first-order scheme, and if there are still  %
% not enough points, the derivative is simply set to zero.                %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'dx'         - REAL SCALAR                                              %
%              - Grid spacing in the 'x' direction.                       %
% ----------------------------------------------------------------------- %
% 'u'          - REAL VECTOR                                              %
%              - One-dimensional real scalar array of data points equally %
%                spaced by 'dx' in the 'x' direction.                     %
% ----------------------------------------------------------------------- %
% 'ux'         - REAL VECTOR                                              %
%              - Derivative of the array u in the 'x' direction.          %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Compute the derivative of the function u = sin(x) on the interval [0,   %
% pi] using a total explicit second-order finite difference scheme with   %
% lower-order switching and a grid spacing of pi/50.                      %
%                                                                         %
% >> dx      = pi/50;                                                     %
% >> u       = sin(0:dx:pi).';                                            %
% >> ux_EXO2 = FD_EXO2(u, dx);                                            %
% >> ux_TRUE = cos(0:dx:pi).';                                            %
% >> E2      = max(abs(ux_EXO2 - ux_TRUE));                               %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% N/A                                                                     %
%                                                                         %
% Called in:                                                              %
% FD_SCH                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FD_EXO2
function [ux] = FD_EXO2(u, dx)

%% Determine the number of data points in the vector u.
n = length(u);

%% If there are more than 2 nodes, use the second-order compact scheme.
ux = zeros(n,1);
if n > 2
    
    %% Compute the derivative at the first boundary node.
    ux(1) = (-u(3) + 4*u(2) - 3*u(1))/(2*dx);
    
    %% Compute the derivative at the last boundary node.
    ux(n) = (u(n-2) - 4*u(n-1) + 3*u(n))/(2*dx);
    
    %% Compute the derivative of the interior nodes.
    for k = 2:n-1
        ux(k) = (u(k+1) - u(k-1))/(2*dx);
    end

%% If there are 2 nodes, switch to a first-order scheme.
elseif n == 2
    
    ux(1) = (u(2) - u(1))/dx;
    ux(2) = ux(1);

%% If there is only one node or none, set the derivative to zero.
else
    
    ux = 0;
    
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
% Line(s) 93                                                              %
% * The scheme being used here is second-order both at the interior nodes %
%   and at the boundary nodes. If only 2 nodes are present, the scheme    %
%   becomes degenerate as the resulting equations become linearly         %
%   dependent.                                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
