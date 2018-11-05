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
% Last Update: October 11th, 2018 by Giuseppe Di Labbio                   %
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
% FTLE = PP_FTLE(VEC, t, dt, rf, freq, dir, options);                     %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the finite-time Lyapunov exponent (FTLE) field   %
% of a discrete velocity field time series evolving with a fixed time     %
% step (or frequency); see [1-4]. While the FTLE field is most commonly   %
% used as a heuristic means of detecting Lagrangian coherent structures   %
% (LCSs), some care must be taken in such a context to consider false     %
% positives and false negatives that may arise [5]. For a review of LCSs  %
% and how the theory pulled away from FTLEs, refer to [6] and the         %
% references therein. This function permits the user to compute the FTLE  %
% field for a specified subset of the time series by integrating forward  %
% or backward over an integration time specified by the user. The user    %
% also has the option of spatially refining the grid using a uniform      %
% refinement factor. Options are made available to use mask information   %
% as well as to append part of the time series ahead or behind the        %
% dataset in the case of periodic flows.                                  %
%                                                                         %
% References:                                                             %
% [1] Haller, G. (2000). Finding finite-time invariant manifolds in two-  %
%     dimensional velocity fields. Chaos, 10(1), 99-108.                  %
% [2] Haller, G. (2001). Distinguished material surfaces and coherent     %
%     structures in three-dimensional fluid flows. Physica D, 149, 248-   %
%     277.                                                                %
% [3] Shadden, S. C., Lekien, F. & Marsden, J. E. (2005). Definition and  %
%     properties of Lagrangian coherent structures from finite-time       %
%     Lyapunov exponents in two-dimensional aperiodic flows. Physica D,   %
%     212, 271-304.                                                       %
% [4] Green, M. A., Rowley, C. W., & Haller, G. (2007). Detection of      %
%     Lagrangian coherent structures in three-dimensional turbulence.     %
%     Journal of Fluid Mechanics, 572, 111-120.                           %
% [5] Haller, G. (2011). A variational theory of hyperbolic Lagrangian    %
%     coherent structures. Physica D, 240, 574-598.                       %
% [6] Haller, G. (2015). Lagrangian coherent structures. Annual Review of %
%     Fluid Mechanics, 47, 137-162.                                       %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'dir'        - STRING                                                   %
%              - The integration direction to compute the finite-time     %
%                Lyapunov exponent field, namely backward ('bwd') or      %
%                forward ('fwd').                                         %
% ----------------------------------------------------------------------- %
% 'dt'         - INTEGER SCALAR                                           %
%              - The number of time steps over which to compute the       %
%                finite-time Lyapunov exponent (the integration time).    %
% ----------------------------------------------------------------------- %
% 'freq'       - REAL SCALAR                                              %
%              - The temporal frequency of the dataset (the reciprocal    %
%                of the duration of each time step).                      %
% ----------------------------------------------------------------------- %
% 'FTLE'       - 1D CELL ARRAY                                            %
%              ---> ELEMENTS: STRUCTS                                     %
%              - Finite-time Lyapunov exponent field stored in a one-     %
%                dimensional cell array of structs. The cell number       %
%                denotes the time step number beginning from the start of %
%                the desired range. Each struct contains as fields the,   %
%                possibly refined, mask ('C'), grid ('X' and 'Y') and     %
%                FTLE field ('Val').                                      %
% ----------------------------------------------------------------------- %
% 'rf'         - INTEGER SCALAR                                           %
%              - The refinement factor denotes the multiple over which    %
%                the grid is to be refined.                               %
% ----------------------------------------------------------------------- %
% 't'          - TWO-ELEMENT INTEGER VECTOR                               %
%              - The range over which to compute the finite-time Lyapunov %
%                exponent in the form [t(1) t(2)].                        %
% ----------------------------------------------------------------------- %
% 'VEC'        - 1D CELL ARRAY                                            %
%              ---> ELEMENTS: STRUCTS                                     %
%              - Velocity field stored in a one-dimensional cell array of %
%                structs. The cell number denotes the time step number.   %
%                Each struct contains as fields the spatial mask ('C'),   %
%                Cartesian grid ('X' and 'Y') and velocity components     %
%                ('U' and 'V').                                           %
% ----------------------------------------------------------------------- %
% Options:                                                                %
% ----------------------------------------------------------------------- %
% 'append'     - SPECIFIC STRING                                          %
%              - Option to append data ahead or behind of the velocity    %
%                field time series by considering the dataset 'VEC' to be %
%                periodic and comprised of one full period. The data is   %
%                appended as many times as needed to satisfy the forward  %
%                or backward integration time defined by 'dt'.            %
%              - Default: Not Active                                      %
% ----------------------------------------------------------------------- %
% 'mask'       - SPECIFIC STRING                                          %
%              - Option to use the mask information contained in the 'C'  %
%                field of the 'VEC' structs.                              %
%              - Default: Not Active                                      %
% ----------------------------------------------------------------------- %
% 'XYinAll'    - SPECIFIC STRING                                          %
%              - Option to store the 'X' and 'Y' fields in all the cells  %
%                of the FTLE variable. If not specified, these fields are %
%                only stored in the first cell.                           %
%              - Default: Not Active                                      %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Calculate the backward finite-time Lyapunov exponent of a time-         %
% dependent double gyre on the domain (x,y) = [0,2]x[0,1] with a constant %
% grid spacing of 0.01 over the time interval [0,20] with time-step size  %
% 0.1. Use A = 0.1, epsilon = 0.25, and omega = 2*pi/10. Use a refinement %
% factor of 8 to calculate the backward FTLE field over the data range    %
% [1 101] with an integration time of 151. Note that you will have to     %
% append data in such a case.                                             %
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

%% Parse the inputs.
if nargin < 6
    error('InputError:tooLittle', 'Not enough input arguments.');
elseif nargin > 5 && nargin < 10
    if ischar(dir)
        chk = {'bwd' 'fwd'}.';
        if ~max(strcmpi(dir, chk))
            error('InputError:notValidDir',                             ...
                 ['The specified integration direction could only be '  ...
                  'one of the following: bwd or fwd.']);
        end
    else
        error('InputError:notString1',                                  ...
              'The integration direction must be a string.');
    end
    tmp = [0 0 0].';
    for k = 1:nargin-6
        chk = {'append' 'mask' 'XYinAll'}.';
        if ischar(varargin{k})
            if ~max(strcmpi(varargin{k}, chk))
                error('InputError:notValidOpt',                             ...
                    ['The specified input options could only be '       ...
                     'selected from the following: append, mask or '    ...
                     'XYinAll.']);
            end
            tmp = max(tmp,strcmpi(varargin{k}, chk));
        else
            error('InputError:notString2', 'The options must be strings.');
        end
    end
    useAppend  = tmp(1);
    useMask    = tmp(2);
    useXYinAll = tmp(3);
    clear chk tmp;
else
    error('InputError:tooMany', 'Too many input arguments.');
end

%% Determine the new grid size in x and y.
nx = rf*size(VEC{1}.X, 2);
ny = rf*size(VEC{1}.Y, 1);

%% Define a refined grid to compute the FTLE field over.
x = linspace(VEC{1}.X(1,1), VEC{1}.X(1,end), nx);
y = linspace(VEC{1}.Y(1,1), VEC{1}.Y(end,1), ny);
[P.X, P.Y] = meshgrid(x, y);
dx = P.X(1,2) - P.X(1,1);
dy = P.Y(1,1) - P.Y(2,1);
clear nx ny x y;

%% Initialize the FTLE data for all times in the time range.
FTLE = cell(t(2)-t(1)+1,1);
if useXYinAll
    for k = 1:length(FTLE)
        FTLE{k} = struct('C', 1, 'Val', 0, 'X', P.X, 'Y', P.Y);
    end
else
    FTLE{1} = struct('C', 1, 'Val', 0, 'X', P.X, 'Y', P.Y);
    for k = 2:length(FTLE)
        FTLE{k} = struct('C', 1, 'Val', 0);
    end
end
clear useXYinAll;

%% Interpolate the mask onto the finer grid.
if useMask
    for k = 1:length(VEC)
        VEC{k}.C = interp2(VEC{k}.X, VEC{k}.Y, double(VEC{k}.C), P.X, P.Y);
        VEC{k}.C = logical(floor(VEC{k}.C));
    end
    for k = 1:length(FTLE)
        FTLE{k}.C = VEC{k+t(1)-1}.C;
    end
end

%% Append data for FTLE calculation.
if useAppend
    if strcmpi(dir, 'bwd')
        VEC = [VEC(end-dt+1:end); VEC];
        t = t + dt;
    elseif strcmpi(dir, 'fwd')
        VEC = [VEC; VEC(1:dt)];
    end
end

if strcmpi(dir, 'bwd') && useMask
    %% Calculation of backward FTLE with a mask.
    
    for k = t(1):t(2)
        %% Advect the particles from k to k-dt.
        A = PP_AdvectGrid(P, flipud(VEC(k-dt:k)), -1/freq, 'singlestep');
        
        %% Compute the FTLE field sequentially over the range k-dt to k.
        for l = 1:dt
            
            % Compute the deformation gradient tensor.
            DG = PP_DeformGradTensor(A{l+1}, [dx dy], 'EXO2', VEC{k-l}.C);
            
            % Compute the elements of the Cauchy-Green strain tensor.
            CG11 = DG.xx.^2 + DG.yx.^2;
            CG12 = DG.xx.*DG.xy + DG.yx.*DG.yy;
            CG22 = DG.xy.^2 + DG.yy.^2;
            
            % Compute the largest eigenvalue of the tensor eigensystem.
            b = CG11 + CG22;
            c = CG11.*CG22 - CG12.^2;
            d = sqrt(b.^2 - 4*c);
            EGV = 0.5*max(b + d, b - d);
            
            % Compute the FTLE from the eigenvalue.
            FTLEnew = log(EGV)/(2*l/freq);
            
            % Replace new NaN values with last known FTLE values.
            if l > 1
                FTLEnew(isnan(FTLEnew)) = FTLEold(isnan(FTLEnew));
            else
                FTLEnew(isnan(FTLEnew)) = 0;
            end
            FTLEold = FTLEnew;
            
        end
        FTLE{k-t(1)+1,1}.Val = FTLEnew;
        
        %% Display the time step count.
        if useAppend
            fprintf('Time Step %d Complete\n', k - dt);
        else
            fprintf('Time Step %d Complete\n', k);
        end
        
    end
elseif strcmpi(dir, 'bwd') && ~useMask
    %% Calculation of backward FTLE without a mask.
    
    for k = t(1):t(2)
        %% Advect the particles from k to k-dt.
        A = PP_AdvectGrid(P, flipud(VEC(k-dt:k)), -1/freq, 'singlestep');
        
        %% Compute the FTLE field sequentially over the range k-dt to k.
        for l = 1:dt
            
            % Compute the deformation gradient tensor.
            DG = PP_DeformGradTensor(A{l+1}, [dx dy], 'EXO2');
            
            % Compute the elements of the Cauchy-Green strain tensor.
            CG11 = DG.xx.^2 + DG.yx.^2;
            CG12 = DG.xx.*DG.xy + DG.yx.*DG.yy;
            CG22 = DG.xy.^2 + DG.yy.^2;
            
            % Compute the largest eigenvalue of the tensor eigensystem.
            b = CG11 + CG22;
            c = CG11.*CG22 - CG12.^2;
            d = sqrt(b.^2 - 4*c);
            EGV = 0.5*max(b + d, b - d);
            
            % Compute the FTLE from the eigenvalue.
            FTLEnew = log(EGV)/(2*l/freq);
            
            % Replace new NaN values with last known FTLE values.
            if l > 1
                FTLEnew(isnan(FTLEnew)) = FTLEold(isnan(FTLEnew));
            else
                FTLEnew(isnan(FTLEnew)) = 0;
            end
            FTLEold = FTLEnew;
            
        end
        FTLE{k-t(1)+1,1}.Val = FTLEnew;

        %% Display the time step count.
        if useAppend
            fprintf('Time Step %d Complete\n', k - dt);
        else
            fprintf('Time Step %d Complete\n', k);
        end
        
    end
elseif strcmpi(dir, 'fwd') && useMask
    %% Calculation of forward FTLE with a mask.
    
    for k = t(1):t(2)
        %% Advect the particles from k to k+dt.
        A = PP_AdvectGrid(P, VEC(k:k+dt), 1/freq, 'singlestep');
        
        %% Compute the FTLE field sequentially over the range k to k+dt.
        for l = 1:dt
            
            % Compute the deformation gradient tensor.
            DG = PP_DeformGradTensor(A{l+1}, [dx dy], 'EXO2', VEC{k+l}.C);
            
            % Compute the elements of the Cauchy-Green strain tensor.
            CG11 = DG.xx.^2 + DG.yx.^2;
            CG12 = DG.xx.*DG.xy + DG.yx.*DG.yy;
            CG22 = DG.xy.^2 + DG.yy.^2;
            
            % Compute the largest eigenvalue of the tensor eigensystem.
            b = CG11 + CG22;
            c = CG11.*CG22 - CG12.^2;
            d = sqrt(b.^2 - 4*c);
            EGV = 0.5*max(b + d, b - d);
            
            % Compute the FTLE from the eigenvalue.
            FTLEnew = log(EGV)/(2*l/freq);
            
            % Replace new NaN values with last known FTLE values.
            if l > 1
                FTLEnew(isnan(FTLEnew)) = FTLEold(isnan(FTLEnew));
            else
                FTLEnew(isnan(FTLEnew)) = 0;
            end
            FTLEold = FTLEnew;
        end
        FTLE{k-t(1)+1,1}.Val = FTLEnew;
        
        %% Display the time step count.
        fprintf('Time Step %d Complete\n', k);
        
    end
elseif strcmpi(dir, 'fwd') && ~useMask
    %% Calculation of forward FTLE without a mask.
    
    for k = t(1):t(2)
        %% Advect the particles from k to k+dt.
        A = PP_AdvectGrid(P, VEC(k:k+dt), 1/freq, 'singlestep');
        
        %% Compute the FTLE field sequentially over the range k to k+dt.
        for l = 1:dt
            
            % Compute the deformation gradient tensor.
            DG = PP_DeformGradTensor(A{l+1}, [dx dy], 'EXO2');
            
            % Compute the elements of the Cauchy-Green strain tensor.
            CG11 = DG.xx.^2 + DG.yx.^2;
            CG12 = DG.xx.*DG.xy + DG.yx.*DG.yy;
            CG22 = DG.xy.^2 + DG.yy.^2;
            
            % Compute the largest eigenvalue of the tensor eigensystem.
            b = CG11 + CG22;
            c = CG11.*CG22 - CG12.^2;
            d = sqrt(b.^2 - 4*c);
            EGV = 0.5*max(b + d, b - d);
            
            % Compute the FTLE from the eigenvalue.
            FTLEnew = log(EGV)/(2*l/freq);
            
            % Replace new NaN values with last known FTLE values.
            if l > 1
                FTLEnew(isnan(FTLEnew)) = FTLEold(isnan(FTLEnew));
            else
                FTLEnew(isnan(FTLEnew)) = 0;
            end
            FTLEold = FTLEnew;
        end
        FTLE{k-t(1)+1,1}.Val = FTLEnew;
        
        %% Display the time step count.
        fprintf('Time Step %d Complete\n', k);
        
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% SUPPRESS MESSAGES %%%%%%%%%%%%%%%%%%%%%%%%% %%

%#ok<*N/A>
% Line(s) N/A
% Message(s)
% * N/A.
% Reason(s)
% * N/A.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Line(s) N/A                                                             %
% * N/A.                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
