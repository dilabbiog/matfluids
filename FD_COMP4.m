%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            FINITE DIFFERENCES                           %
%                     LELE'S 4TH ORDER COMPACT SCHEME                     %
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
% ux = FD_COMP4(u, dx);                                                   %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function applies the classical fourth-order Padé (or compact)      %
% finite difference scheme described in [1] for the first derivative. The %
% boundary node scheme is taken as a fourth-order scheme as well, which   %
% is also given in [1]. In cases where there are not enough points to     %
% apply the total fourth-order scheme, this function switches to lower-   %
% order schemes.                                                          %
%                                                                         %
% References:                                                             %
% [1] Lele, S. K. (1992). Compact finite difference schemes with          %
%     spectral-like resolution. Journal of Computational Physics, 103,    %
%     16-42.                                                              %
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
% pi] using a total compact fourth-order finite difference scheme with    %
% lower-order switching and a grid spacing of pi/50.                      %
%                                                                         %
% >> dx      = pi/50;                                                     %
% >> u       = sin(0:dx:pi).';                                            %
% >> ux_COM4 = FD_COMP4(u, dx);                                           %
% >> ux_TRUE = cos(0:dx:pi).';                                            %
% >> E4      = max(abs(ux_COM4 - ux_TRUE));                               %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% FD_COMP2                                                                %
% FD_COMP3                                                                %
% TDMA                                                                    %
%                                                                         %
% Called in:                                                              %
% FD_SCH                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FD_COMP4
function [ux] = FD_COMP4(u, dx)

%% Determine the number of data points in the vector u.
n = length(u);

%% If there are more than 4 nodes, use the fourth-order compact scheme.
if n > 4
    
    %% Define the scheme constants for the interior nodes.
    alpha_i = 1/4;
    a_i = 3/2;
    
    %% Define the the scheme constants for boundary nodes.
    alpha_b = 3;
    a_b = -17/6;
    b_b = 3/2;
    c_b = 3/2;
    d_b = -1/6;
    
    %% Define the lower and upper diagonals for the tridiagonal system.
    LD = [alpha_i*ones(n-2,1); alpha_b];
    UD = [alpha_b; alpha_i*ones(n-2,1)];
    
    %% Initialize the column vector b from [A][ux] = [b].
    b = zeros(n,1);
    
    %% Define the entries of the column vector b.
    for k = 1:n
        if k == 1
            % If the node is the first boundary node, use the appropriate
            % fourth-order one-way scheme.
            b(k) = (a_b*u(k) + b_b*u(k+1) + c_b*u(k+2) + d_b*u(k+3))/dx;
        elseif k == n
            % If the node is the last boundary node, use the appropriate
            % fourth-order one-way scheme.
            b(k) = -(a_b*u(k) + b_b*u(k-1) + c_b*u(k-2) + d_b*u(k-3))/dx;
        else
            % If the node is an interior node, use the appropriate fourth-
            % order central scheme.
            b(k) = a_i*(u(k+1) - u(k-1))/(2*dx);
        end
    end
    
    %% Compute the derivative using the tridiagonal matrix algorithm.
    ux = TDMA(LD, ones(n,1), UD, b);

%% If there are 4 nodes or less, switch to the third-order compact scheme.
else
    
    ux = FD_COMP3(u, dx);
    
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
% Line(s) 96                                                              %
% * The scheme being used here is fourth-order both at the interior nodes %
%   and at the boundary nodes. If only 4 nodes are present, the scheme    %
%   becomes degenerate as the resulting equations become linearly         %
%   dependent.                                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
