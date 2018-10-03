%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            FINITE DIFFERENCES                           %
%                SCHEMES USED WITH OR WITHOUT GEOMETRY DATA               %
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
% ux = FD_SCH(u, dx, scheme);                                             %
% ux = FD_SCH(u, dx, scheme, geo);                                        %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the derivative of a one-dimensional array by     %
% finite-differences. The user has the option of selecting which finite   %
% difference scheme will be used.                                         %
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
% 'geo'        - REAL VECTOR                                              %
%              - One-dimensional (scalar) array of points containing mask %
%                information for the array u. Values of 1 correspond to   %
%                active nodes while values of 0 correspond to blanked     %
%                nodes.                                                   %
% ----------------------------------------------------------------------- %
% 'scheme'     - SPECIFIC STRING                                          %
%              - String denoting which scheme to use. The user can select %
%                from the explicit second-order central scheme ('EXO2'),  %
%                or Lele's second, third or fourth order compact schemes  %
%                ('COMP2', 'COMP3', 'COMP4').                             %
% ----------------------------------------------------------------------- %
% 'u'          - REAL VECTOR                                              %
%              - One-dimensional real scalar array of data points equally %
%                spaced by 'dx' in the 'x' direction.                     %
% ----------------------------------------------------------------------- %
% 'ux'         - REAL VECTOR                                              %
%              - Derivative of the array u in the 'x' direction.          %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Compute the derivative of the function u = sin(x) on the interval [0,   %
% pi] using a total explicit second-order finite difference scheme with   %
% lower-order switching and a grid spacing of pi/50.                      %
%                                                                         %
% >> dx      = pi/50;                                                     %
% >> u       = sin(0:dx:pi).';                                            %
% >> ux_EXO2 = FD_SCH(u, dx, 'EXO2');                                     %
% >> ux_TRUE = cos(0:dx:pi).';                                            %
% >> E2      = max(abs(ux_EXO2 - ux_TRUE));                               %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Compute the derivative of the function u = sin(x) on the interval [0,   %
% pi] using a total explicit second-order finite difference scheme with   %
% lower-order switching and a grid spacing of pi/50 for some mask g.      %
%                                                                         %
% >> dx  = pi/50;                                                         %
% >> geo = [ones(20,1); zeros(10,1); ones(6,1); zeros(4,1); ones(11,1)];  %
% >> u   = geo.*(sin(0:dx:pi).');                                         %
% >> ux_EXO2 = FD_SCH(u, dx, 'EXO2', geo);                                %
% >> ux_TRUE = geo.*(cos(0:dx:pi).');                                     %
% >> E2      = max(abs(ux_EXO2 - ux_TRUE));                               %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% FD_COMP2                                                                %
% FD_COMP3                                                                %
% FD_COMP4                                                                %
% FD_EXO2                                                                 %
% TDMA                                                                    %
%                                                                         %
% Called in:                                                              %
% FD_CR4                                                                  %
% PP_CauchyGreenTensor                                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FD_SCH
function [ux] = FD_SCH(u, dx, scheme, varargin)

if nargin == 3
    % Depending on the option, apply the appropriate finite difference
    % scheme.
    
    if strcmpi(scheme,'COMP2')
        ux = FD_COMP2(u, dx);
    elseif strcmpi(scheme,'COMP3')
        ux = FD_COMP3(u, dx);
    elseif strcmpi(scheme,'COMP4')
        ux = FD_COMP4(u, dx);
    elseif strcmpi(scheme,'EXO2')
        ux = FD_EXO2(u, dx);
    end
    
elseif nargin == 4
    % Determine the number of nodes (n) in the vector (u).
    n = length(u);
    
    % Initialize the derivative vector (ux) of the vector (u).
    ux = zeros(n,1);
    
    % Initialize a counter to determine sets of nodes by isolation.
    count = 0;
    
    for k = 1:n
        
        if ~varargin{1}(k)
            % If the current node is an exterior node (= 0), continue to
            % the next node.
            continue;
        else
            % If the current node is a non-exterior node, increment the
            % counter for the number of isolated nodes.
            count = count + 1;
            
            % If the end of the domain is reached (and is a non-exterior
            % node) or if the next node is an exterior node, compute the
            % derivative vector using the appropriate finite-difference
            % scheme.
            if k == n || ~varargin{1}(k+1)
                % Depending on the option, apply the appropriate finite
                % difference scheme.
                if strcmpi(scheme,'COMP2')
                    ux(k-count+1:k) = FD_COMP2(u(k-count+1:k), dx);
                elseif strcmpi(scheme,'COMP3')
                    ux(k-count+1:k) = FD_COMP3(u(k-count+1:k), dx);
                elseif strcmpi(scheme,'COMP4')
                    ux(k-count+1:k) = FD_COMP4(u(k-count+1:k), dx);
                elseif strcmpi(scheme,'EXO2')
                    ux(k-count+1:k) = FD_EXO2(u(k-count+1:k), dx);
                end
                
                % Reset the isolated node counter.
                count = 0;
            end
        end
    end
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
% Line(s) N/A                                                             %
% * N/A.                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
