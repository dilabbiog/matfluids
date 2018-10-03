%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                   INTERPOLATE VELOCITY FIELD IN SPACE                   %
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
% iVEC = PP_InterpVelFieldSpace(VEC, N, method);                          %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function interpolates the vector field in a refined spatial domain %
% defined by some refinement factor N.                                    %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'iVEC'       - STRUCT                                                   %
%              - Refined velocity vector field in the form of a struct    %
%                containing information on the spatial mask (C), velocity %
%                components (U and V), and Cartesian grid (X and Y).      %
% ----------------------------------------------------------------------- %
% 'method'     - SPECIFIC STRING                                          %
%              - Method of interpolation to use. See MATLAB's interp2     %
%                function for a list of options.                          %
% ----------------------------------------------------------------------- %
% 'N'          - INTEGER SCALAR                                           %
%              - Refinement factor of the spatial field.                  %
% ----------------------------------------------------------------------- %
% 'VEC'        - STRUCT                                                   %
%              - Struct containing information on the spatial mask (C),   %
%                velocity components (U and V), and Cartesian grid (X and %
%                Y).                                                      %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Refine the time series of a time-dependent double gyre in space after   %
% creating one over the domain (x,y) = [0,2]x[0,1] with a constant grid   %
% spacing of 0.01 over the time interval [0,20] with time-step size 0.1.  %
% Use A = 0.1, epsilon = 0.25, and omega = 2*pi/10. Use a refinement      %
% factor of 4 and use the 'cubic' method (i.e. cubic convolution).        %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A    = 0.1;                                                          %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> iVEC = cell(length(VEC),1);                                          %
% >> for k = 1:length(t)                                                  %
% iVEC{k} = PP_InterpVelFieldSpace(VEC{k}, 4, 'cubic');                   %
% end                                                                     %
% >> for k = 1:length(t)                                                  %
% quiver(iVEC{k}.X(1:16:end,1:16:end), iVEC{k}.Y(1:16:end,1:16:end), ...  %
%        iVEC{k}.U(1:16:end,1:16:end), iVEC{k}.V(1:16:end,1:16:end));     %
% axis([0 2 0 1]);                                                        %
% pause(0.25);                                                            %
% end                                                                     %
% >> clear k;                                                             %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% N/A                                                                     %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_InterpVelFieldSpace
function [iVEC] = PP_InterpVelFieldSpace(VEC, N, method)

% Initialize the interpolated vector field struct.
iVEC = VEC;

% Determine the new grid size in x and y.
nx = N*size(VEC.X, 2);
ny = N*size(VEC.Y, 1);

% Define a refined grid for the interpolation.
x = linspace(VEC.X(1,1), VEC.X(1,end), nx);
y = fliplr(linspace(VEC.Y(end,1), VEC.Y(1,1), ny));
[iVEC.X, iVEC.Y] = meshgrid(x, y);
clear nx ny x y;

% Interpolate the U and V components of velocity.
iVEC.U = interp2(VEC.X, VEC.Y, VEC.U, iVEC.X, iVEC.Y, method);
iVEC.V = interp2(VEC.X, VEC.Y, VEC.V, iVEC.X, iVEC.Y, method);

% Interpolate for the new mask.
iVEC.C = interp2(VEC.X, VEC.Y, double(VEC.C), iVEC.X, iVEC.Y, 'linear');
iVEC.C = logical(iVEC.C);

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
