%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                     LOCAL ACCELERATION VECTOR FIELD                     %
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
% ACC = PP_LocalAccel(VEC, dt);                                           %
% ACC = PP_LocalAccel(VEC, dt, 'mask');                                   %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the local acceleration vector field of a         %
% velocity field. The function detects if the calculation is two or three %
% dimensional.                                                            %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'ACC'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - Local acceleration vector field time series in the form  %
%                of a cell array of structs. The cell number denotes the  %
%                time step number. Each cell is a struct containing the   %
%                vector components of the acceleration field ('x' and 'y' %
%                in two dimensions). In three dimensions, the structs     %
%                ought to contain the components 'x', 'y' and 'z'.        %
% ----------------------------------------------------------------------- %
% 'dt'         - REAL SCALAR                                              %
%              - Time step of the raw data set.                           %
% ----------------------------------------------------------------------- %
% 'opt'        - SPECIFIC STRING                                          %
%              - Option to use a pre-defined mask in VEC.C using 'mask'   %
%                to both speed up the calculation and treat the boundary  %
%                nodes correctly.                                         %
% ----------------------------------------------------------------------- %
% 'VEC'        - 1D CELL, ELEMENTS: STRUCTS (C, X, Y, U, V)               %
%              - Velocity vector field time series in the form of a cell  %
%                array of structs. The cell number denotes the time step  %
%                number. Each cell is a struct containing the mask ('C'), %
%                Cartesian coordinates ('X' and 'Y') and velocity vector  %
%                components ('U' and 'V'). In three dimensions each cell  %
%                struct contains the mask ('C'), Cartesian coordinates    %
%                ('X', 'Y' and 'Z') and velocity vector components ('U',  %
%                'V' and 'W').                                            %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Calculate the local acceleration vector field of a time-dependent       %
% double gyre on the domain (x,y) = [0,2]x[0,1] with a constant grid      %
% spacing of 0.01 over the time interval [0,20] with time-step size 0.1.  %
% Use A = 0.1, epsilon = 0.25, and omega = 2*pi/10.                       %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> ACC = PP_LocalAccel(VEC, t(2));                                      %
% >> Amag = cell(length(t),1);                                            %
% >> for k = 1:length(t)                                                  %
% Amag{k} = sqrt(ACC{k}.x.^2 + ACC{k}.y.^2);                              %
% end                                                                     %
% >> for k = 1:length(t)                                                  %
% contourf(VEC{k}.X, VEC{k}.Y, Amag{k}, 'EdgeColor', 'None');             %
% colormap jet;                                                           %
% hold on;                                                                %
% quiver(VEC{k}.X(1:4:end,1:4:end), VEC{k}.Y(1:4:end,1:4:end), ...        %
%        VEC{k}.U(1:4:end,1:4:end), VEC{k}.V(1:4:end,1:4:end), 'k');      %
% hold off;                                                               %
% pause(0.25);                                                            %
% end                                                                     %
% >> clear k;                                                             %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% FD_COMP2                                                                %
% FD_COMP3                                                                %
% FD_COMP4                                                                %
% FD_EXO2                                                                 %
% FD_SCH                                                                  %
% field2mat3D                                                             %
% TDMA                                                                    %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_LocalAccel
function [ACC] = PP_LocalAccel(VEC, dt, varargin)

if length(fieldnames(VEC{1})) == 5
    
    [r,c] = size(VEC{1}.X);
    
    ACC = cell(length(VEC),1);
    for k = 1:length(VEC)
        ACC{k} = struct('x', 0, 'y', 0);
    end
    
    U = field2mat3D(VEC, 'U');
    V = field2mat3D(VEC, 'V');
    
    Ut = zeros(r,c,length(VEC));
    Vt = zeros(r,c,length(VEC));
    if nargin == 2
        for i = 1:r
            for j = 1:c
                Ut(i,j,:) = FD_SCH(U(i,j,:), dt, 'COMP4');
                Vt(i,j,:) = FD_SCH(V(i,j,:), dt, 'COMP4');
            end
        end
    elseif nargin == 3 && strcmpi(varargin{1},'mask')
        C = field2mat3D(VEC, 'C');
        for i = 1:r
            for j = 1:c
                if ~max(C(i,j,:))
                    continue;
                else
                    Ut(i,j,:) = FD_SCH(U(i,j,:), dt, 'COMP4');
                    Vt(i,j,:) = FD_SCH(V(i,j,:), dt, 'COMP4');
                end
            end
        end
    end
    
    for k = 1:length(VEC)
        ACC{k}.x = Ut(:,:,k);
        ACC{k}.y = Vt(:,:,k);
    end
    
elseif length(fieldnames(VEC{1})) == 7
    
    % Coming soon...
    
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
