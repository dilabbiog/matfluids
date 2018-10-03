%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                      FINITE-TIME LYAPUNOV EXPONENT                      %
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
% FTLE = PP_FTLE(VEC, t, dt, rf, freq, dir);                              %
% FTLE = PP_FTLE(VEC, t, dt, rf, freq, dir, 'append');                    %
% FTLE = PP_FTLE(VEC, t, dt, rf, freq, dir, 'mask');                      %
% FTLE = PP_FTLE(VEC, t, dt, rf, freq, dir, 'mask', 'append');            %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the forward or backward finite-time Lyapunov     %
% exponent of vector field evolving in time (with a known frequency). The %
% range, step range and refinement factor are all adjustable. The user    %
% has the option of using mask information as well as of appending part   %
% of the time series ahead or behind the dataset in the case of periodic  %
% flows.                                                                  %
%                                                                         %
% References:                                                             %
% Coming soon ...                                                         %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'dir'        - SPECIFIC STRING                                          %
%              - The direction in which to compute the finite-time        %
%                Lyapunov exponent, namely backward ('bwd') or forward    %
%                ('fwd').                                                 %
% ----------------------------------------------------------------------- %
% 'dt'         - INTEGER SCALAR                                           %
%              - The number of time steps over which to compute the       %
%                finite-time Lyapunov exponent.                           %
% ----------------------------------------------------------------------- %
% 'freq'       - REAL SCALAR                                              %
%              - The temporal frequency of the data set (the reciprocal   %
%                of the duration of each time step).                      %
% ----------------------------------------------------------------------- %
% 'FTLE'       - 1D CELL, ELEMENTS: STRUCTS                               %
%              - Finite-time Lyapunov exponent field stored in a cell     %
%                array of structs. The cell number denotes the time step  %
%                number beginning from the desired range. Each cell is a  %
%                struct containing the refined mask ('C'), the refined    %
%                grid ('X' and 'Y') and FTLE value ('Val').               %
% ----------------------------------------------------------------------- %
% 'opt'        - SPECIFIC STRING                                          %
%              - Option to use the mask information contained in VEC.C    %
%                using 'mask' and to append data ahead or behind of the   %
%                velocity vector field using 'append'.                    %
% ----------------------------------------------------------------------- %
% 'rf'         - INTEGER SCALAR                                           %
%              - The refinement factor denotes the multiple over which    %
%                the grid is to be refined.                               %
% ----------------------------------------------------------------------- %
% 't'          - TWO-ELEMENT INTEGER VECTOR                               %
%              - The range over which to compute the finite-time Lyapunov %
%                exponent in the form [t(1) t(2)].                        %
% ----------------------------------------------------------------------- %
% 'VEC'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - One-dimensional cell array of structs each containing    %
%                information on the spatial mask (C), velocity components %
%                (U and V), and Cartesian grid (X and Y).                 %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Calculate the backward finite-time Lyapunov exponent of a time-         %
% dependent double gyre on the domain (x,y) = [0,2]x[0,1] with a constant %
% grid spacing of 0.01 over the time interval [0,20] with time-step size  %
% 0.1. Use A = 0.1, epsilon = 0.25, and omega = 2*pi/10. Use a refinement %
% factor of 8 to calculate the backward FTLE field over the full data     %
% range (i.e. [1 21]) with an integration step of 10. Note that you will  %
% have to append data in such a case.                                     %
%                                                                         %
% >> x = linspace(0, 2, 41).';                                            %
% >> y = linspace(0, 1, 21).';                                            %
% >> t = linspace(0, 20, 201).';                                          %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> FTLE = PP_FTLE(VEC, [1 101], 151, 8, 1/t(2), 'fwd', 'append');       %
% >> for k = 1:length(FTLE)                                               %
% contourf(FTLE{k}.X, FTLE{k}.Y, FTLE{k}.Val, 'EdgeColor', 'None');       %
% colormap gray;                                                          %
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
% PP_AdvectGrid                                                           %
% PP_DeformGradTensor                                                     %
% TDMA                                                                    %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_FTLE
function [FTLE] = PP_FTLE(VEC, t, dt, rf, freq, dir, varargin)

% Determine the new grid size in x and y.
nx = rf*size(VEC{1}.X, 2);
ny = rf*size(VEC{1}.Y, 1);

% Define a refined grid to compute the FTLE field over.
x = linspace(VEC{1}.X(1,1), VEC{1}.X(1,end), nx);
y = fliplr(linspace(VEC{1}.Y(end,1), VEC{1}.Y(1,1), ny));
[P.X, P.Y] = meshgrid(x, y);
dx = P.X(1,2) - P.X(1,1);
dy = P.Y(1,1) - P.Y(2,1);
clear x y;

% Initialize the FTLE data for all times in the time range.
FTLE = cell(t(2)-t(1)+1,1);
for k = 1:length(FTLE)
    FTLE{k} = struct('C', 1, 'Val', zeros(ny,nx), 'X', P.X, 'Y', P.Y);
end

% Interpolate the mask onto the finer grid.
flagM = 0;
for v = 1:nargin-6
    if strcmpi(varargin{v}, 'mask')
        flagM = 1;
        for k = 1:length(VEC)
            VEC{k}.C = interp2(VEC{k}.X, VEC{k}.Y, double(VEC{k}.C), ...
                               P.X, P.Y);
            VEC{k}.C = logical(VEC{k}.C);
        end
        for k = 1:length(FTLE)
            FTLE{k}.C = VEC{k+t(1)-1}.C;
        end
        break
    end
end

% Append data for FTLE calculation.
flagA = 0;
for v = 1:nargin-6
    if strcmpi(varargin{v}, 'append')
        flagA = 1;
        if strcmpi(dir, 'bwd')
            VEC = [VEC(end-dt+1:end); VEC];
            t = t + dt;
        elseif strcmpi(dir, 'fwd')
            VEC = [VEC; VEC(1:dt)];
        end
        break
    end
end

if strcmpi(dir, 'bwd') && flagM
    
    % Calculation of backward finite-time Lyapunov exponent with a mask.
    for k = t(1):t(2)
        % Advect the particles from k to k-dt.
        A = PP_AdvectGrid(P, flipud(VEC(k-dt:k)), -1/freq, 'singlestep');
        
        % Compute the FTLE field sequentially over the range k-dt to k.
        CG = cell(dt,1);
        subFTLE = cell(dt,1);
        for l = 1:dt
            CG{l} = PP_DeformGradTensor(A{l+1}, [dx dy], 'EXO2', ...
                                        VEC{k-l}.C);
            EGV2 = ones(ny, nx);
            for i = 1:ny
                for j = 1:nx
                    if ~VEC{k-l}.C(i,j)
                        continue;
                    else
                        M = [CG{l}.xx(i,j) CG{l}.xy(i,j);
                             CG{l}.yx(i,j) CG{l}.yy(i,j)];
                        if ~max(max(isnan(M)))
                            EGV2(i,j) = max(eig(M'*M));
                        else
                            EGV2(i,j) = NaN;
                        end
                    end
                    
                end
            end
            subFTLE{l} = log(EGV2)/(2*l/freq);
            if l > 1
                subFTLE{l}(isnan(subFTLE{l})) = ...
                    subFTLE{l-1}(isnan(subFTLE{l}));
            else
                subFTLE{l}(isnan(subFTLE{l})) = 0;
            end
        end
        FTLE{k-t(1)+1,1}.Val = subFTLE{end};
        
        if flagA
            fprintf('Time Step %d Complete\n', k - dt);
        else
            fprintf('Time Step %d Complete\n', k);
        end
        
    end
elseif strcmpi(dir, 'bwd') && ~flagM
    % Calculation of backward finite-time Lyapunov exponent without a mask.
    for k = t(1):t(2)
        % Advect the particles from k to k-dt.
        A = PP_AdvectGrid(P, flipud(VEC(k-dt:k)), -1/freq, 'singlestep');
        
        % Compute the FTLE field sequentially over the range k-dt to k.
        CG = cell(dt,1);
        subFTLE = cell(dt,1);
        for l = 1:dt
            CG{l} = PP_DeformGradTensor(A{l+1}, [dx dy], 'EXO2');
            EGV2 = ones(ny, nx);
            for i = 1:ny
                for j = 1:nx
                    M = [CG{l}.xx(i,j) CG{l}.xy(i,j);
                         CG{l}.yx(i,j) CG{l}.yy(i,j)];
                    if ~max(max(isnan(M)))
                        EGV2(i,j) = max(eig(M'*M));
                    else
                        EGV2(i,j) = NaN;
                    end
                end
            end
            subFTLE{l} = log(EGV2)/(2*l/freq);
            if l > 1
                subFTLE{l}(isnan(subFTLE{l})) = ...
                    subFTLE{l-1}(isnan(subFTLE{l}));
            else
                subFTLE{l}(isnan(subFTLE{l})) = 0;
            end
        end
        FTLE{k-t(1)+1,1}.Val = subFTLE{end};
        
        if flagA
            fprintf('Time Step %d Complete\n', k - dt);
        else
            fprintf('Time Step %d Complete\n', k);
        end
        
    end
elseif strcmpi(dir, 'fwd') && flagM
    
    % Calculation of forward finite-time Lyapunov exponent with a mask.
    for k = t(1):t(2)
        % Advect the particles from k to k+dt.
        A = PP_AdvectGrid(P, VEC(k:k+dt), 1/freq, 'singlestep');
        
        % Compute the FTLE field sequentially over the range k to k+dt.
        CG = cell(dt,1);
        subFTLE = cell(dt,1);
        for l = 1:dt
            CG{l} = PP_DeformGradTensor(A{l+1}, [dx dy], 'EXO2', ...
                                        VEC{k+l}.C);
            EGV2 = ones(ny, nx);
            for i = 1:ny
                for j = 1:nx
                    if ~VEC{k+l}.C(i,j)
                        continue;
                    else
                        M = [CG{l}.xx(i,j) CG{l}.xy(i,j);
                             CG{l}.yx(i,j) CG{l}.yy(i,j)];
                        if ~max(max(isnan(M)))
                            EGV2(i,j) = max(eig(M'*M));
                        else
                            EGV2(i,j) = NaN;
                        end
                    end
                    
                end
            end
            subFTLE{l} = log(EGV2)/(2*l/freq);
            if l > 1
                subFTLE{l}(isnan(subFTLE{l})) = ...
                    subFTLE{l-1}(isnan(subFTLE{l}));
            else
                subFTLE{l}(isnan(subFTLE{l})) = 0;
            end
        end
        FTLE{k-t(1)+1,1}.Val = subFTLE{end};
        
        fprintf('Time Step %d Complete\n', k);
        
    end
elseif strcmpi(dir, 'fwd') && ~flagM
    
    % Calculation of forward finite-time Lyapunov exponent without a mask.
    for k = t(1):t(2)
        % Advect the particles from k to k+dt.
        A = PP_AdvectGrid(P, VEC(k:k+dt), 1/freq, 'singlestep');
        
        % Compute the FTLE field sequentially over the range k to k+dt.
        CG = cell(dt,1);
        subFTLE = cell(dt,1);
        for l = 1:dt
            CG{l} = PP_DeformGradTensor(A{l+1}, [dx dy], 'EXO2');
            EGV2 = ones(ny, nx);
            for i = 1:ny
                for j = 1:nx
                    M = [CG{l}.xx(i,j) CG{l}.xy(i,j);
                         CG{l}.yx(i,j) CG{l}.yy(i,j)];
                    if ~max(max(isnan(M)))
                        EGV2(i,j) = max(eig(M'*M));
                    else
                        EGV2(i,j) = NaN;
                    end
                end
            end
            subFTLE{l} = log(EGV2)/(2*l/freq);
            if l > 1
                subFTLE{l}(isnan(subFTLE{l})) = ...
                    subFTLE{l-1}(isnan(subFTLE{l}));
            else
                subFTLE{l}(isnan(subFTLE{l})) = 0;
            end
        end
        FTLE{k-t(1)+1,1}.Val = subFTLE{end};
        
        fprintf('Time Step %d Complete\n', k);
        
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% SUPPRESS MESSAGES %%%%%%%%%%%%%%%%%%%%%%%%% %%

%#ok<*AGROW>
% Line(s) N/A
% Message(s)
% * The variable 'VEC' appears to change size on every loop iteration.
%   Consider preallocating for speed.
% Reason(s)
% * This is false since the commands are only executed once.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Line(s) N/A                                                             %
% * N/A.                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
