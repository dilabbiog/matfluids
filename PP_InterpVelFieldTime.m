%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                    INTERPOLATE VELOCITY FIELD IN TIME                   %
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
% iVEC = PP_InterpVelFieldTime(VEC, T, N, method);                        %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function interpolates the vector field in a refined temporal       %
% domain defined by some refinement factor N. Note that the final vector  %
% field time series will have one less than N times the original number   %
% of time steps.                                                          %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'iVEC'       - 1D CELL, ELEMENTS: STRUCTS                               %
%              - Refined velocity vector field in the form of a one-      %
%                dimensional cell array of structs each containing        %
%                information on the spatial mask (C), velocity components %
%                (U and V), and Cartesian grid (X and Y).                 %
% ----------------------------------------------------------------------- %
% 'method'     - SPECIFIC STRING                                          %
%              - Method of interpolation to use. See MATLAB's interp1     %
%                function for a list of options.                          %
% ----------------------------------------------------------------------- %
% 'N'          - INTEGER SCALAR                                           %
%              - Refinement factor of the spatial field.                  %
% ----------------------------------------------------------------------- %
% 'T'          - REAL VECTOR                                              %
%              - Time vector of the data.                                 %
% ----------------------------------------------------------------------- %
% 'VEC'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - One-dimensional cell array of structs each containing    %
%                information on the spatial mask (C), velocity components %
%                (U and V), and Cartesian grid (X and Y).                 %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Refine the time series of a time-dependent double gyre after creating   %
% one over the domain (x,y) = [0,2]x[0,1] with a constant grid spacing of %
% 0.01 over the time interval [0,20] with time-step size 0.1. Use A =     %
% 0.1, epsilon = 0.25, and omega = 2*pi/10. Use a refinement factor of 4  %
% and use the 'pchip' method (i.e. shape-preserving piecewise cubic       %
% interpolation).                                                         %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A    = 0.1;                                                          %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> iVEC = PP_InterpVelFieldTime(VEC, t, 4, 'pchip');                    %
% >> t2 = (0:t(2)/4:t(end)).';                                            %
% >> for k = 1:length(t2)                                                 %
% quiver(iVEC{k}.X(1:4:end,1:4:end), iVEC{k}.Y(1:4:end,1:4:end), ...      %
%        iVEC{k}.U(1:4:end,1:4:end), iVEC{k}.V(1:4:end,1:4:end));         %
% axis([0 2 0 1]);                                                        %
% pause(0.05);                                                            %
% end                                                                     %
% >> clear k;                                                             %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% field2mat3D                                                             %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_InterpVelFieldTime
function [iVEC] = PP_InterpVelFieldTime(VEC, T, N, method)

dt = (T(2) - T(1))/N;
C = permute(field2mat3D(VEC, 'C'), [3 1 2]);
X = permute(field2mat3D(VEC, 'X'), [3 1 2]);
Y = permute(field2mat3D(VEC, 'Y'), [3 1 2]);
U = permute(field2mat3D(VEC, 'U'), [3 1 2]);
V = permute(field2mat3D(VEC, 'V'), [3 1 2]);

Cfrac = permute(interp1(T, double(C), T(1):dt:T(end), method), [2 3 1]);
Xfrac = permute(interp1(T, X, T(1):dt:T(end), method), [2 3 1]);
Yfrac = permute(interp1(T, Y, T(1):dt:T(end), method), [2 3 1]);
Ufrac = permute(interp1(T, U, T(1):dt:T(end), method), [2 3 1]);
Vfrac = permute(interp1(T, V, T(1):dt:T(end), method), [2 3 1]);
clear dt X Y U V;

iVEC = cell(size(Xfrac,3),1);
for k = 1:length(iVEC)
    iVEC{k}.C = logical(Cfrac(:,:,k));
    iVEC{k}.X = Xfrac(:,:,k);
    iVEC{k}.Y = Yfrac(:,:,k);
    iVEC{k}.U = iVEC{k}.C.*Ufrac(:,:,k);
    iVEC{k}.V = iVEC{k}.C.*Vfrac(:,:,k);
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
