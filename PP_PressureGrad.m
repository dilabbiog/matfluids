%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                      PRESSURE GRADIENT VECTOR FIELD                     %
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
% SYNTAX                                                                  %
%                                                                         %
% PG = PP_PressureGrad(VEC, VGT, dt, rho, mu, BF);                        %
% PG = PP_PressureGrad(VEC, VGT, dt, rho, mu, BF, 'mask');                %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the pressure gradient vector field of a velocity %
% field given its velocity gradient tensor, fluid properties and body     %
% force. The function detects if the calculation is two or three          %
% dimensional.                                                            %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'dt'         - REAL SCALAR                                              %
%              - Time step of the raw data set.                           %
% ----------------------------------------------------------------------- %
% 'mu'         - REAL SCALAR                                              %
%              - Dynamic viscosity of the fluid.                          %
% ----------------------------------------------------------------------- %
% 'opt'        - SPECIFIC STRING                                          %
%              - Option to use a pre-defined mask in VEC.C using 'mask'   %
%                to both speed up the calculation and treat the boundary  %
%                nodes correctly.                                         %
% ----------------------------------------------------------------------- %
% 'PG'         - 1D CELL, ELEMENTS: STRUCTS (x, y, z)                     %
%              - Pressure gradient vector field time series in the form   %
%                of a cell array of structs. The cell number denotes the  %
%                time step number. Each cell is a struct containing the   %
%                vector components of the pressure gradient ('x' and 'y'  %
%                in two dimensions). In three dimensions, the structs     %
%                ought to contain the components 'x', 'y' and 'z'.        %
% ----------------------------------------------------------------------- %
% 'rho'        - REAL SCALAR                                              %
%              - Density of the fluid.                                    %
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
% 'VGT'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - Velocity gradient tensor time series in the form of a    %
%                cell array of structs. The cell number denotes the time  %
%                step number. Each cell is a struct containing the four   %
%                components 'UX', 'UY', 'VX' and 'VY' in 2D or the nine   %
%                components 'UX', 'UY', 'UZ', 'VX', 'VY', 'VZ', 'WX',     %
%                'WY' and 'WZ' in 3D.                                     %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Calculate the pressure gradient vector field of a time-dependent double %
% gyre on the domain (x,y) = [0,2]x[0,1] with a constant grid spacing of  %
% 0.01 over the time interval [0,20] with time-step size 0.1. Use A =     %
% 0.1, epsilon = 0.25, and omega = 2*pi/10. Use the density and viscosity %
% of water (1000 kg/m^3 and 0.001 Pa*s respectively). Note that there is  %
% no body force for this flow.                                            %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> mu   = 0.001;                                                        %
% >> rho  = 1000;                                                         %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> BF = cell(length(t),1);                                              %
% >> for k = 1:length(t)                                                  %
% BF{k} = struct('x', 0, 'y', 0);                                         %
% end                                                                     %
% >> PG = PP_PressureGrad(VEC, VGT, t(2), rho, mu, BF);                   %
% >> PGmag = cell(length(t),1);                                           %
% >> for k = 1:length(t)                                                  %
% PGmag{k} = sqrt(PG{k}.x.^2 + PG{k}.y.^2);                               %
% end                                                                     %
% >> for k = 1:length(t)                                                  %
% contourf(VEC{k}.X, VEC{k}.Y, PGmag{k}, 'EdgeColor', 'None');            %
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
% FD_CR4                                                                  %
% FD_EXO2                                                                 %
% FD_SCH                                                                  %
% field2mat3D                                                             %
% PP_AdvectiveAccel                                                       %
% PP_LocalAccel                                                           %
% PP_TotalAccel                                                           %
% TDMA                                                                    %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_PressureGrad
function [PG] = PP_PressureGrad(VEC, VGT, dt, rho, mu, BF, varargin)

if length(fieldnames(VEC{1})) == 5
    
    if nargin == 6
        
        ACC = PP_TotalAccel(VEC, VGT, dt);
        
        UXX = cell(length(VEC),1);
        UYY = cell(length(VEC),1);
        VXX = cell(length(VEC),1);
        VYY = cell(length(VEC),1);
        for k = 1:length(VEC)
            
            TMP1 = struct('C', VEC{k}.C, 'X', VEC{k}.X, 'Y', VEC{k}.Y, ...
                          'U', VGT{k}.UX, 'V', VGT{k}.UY);
            TMP1 = FD_CR4(TMP1);
            TMP2 = struct('C', VEC{k}.C, 'X', VEC{k}.X, 'Y', VEC{k}.Y, ...
                          'U', VGT{k}.VX, 'V', VGT{k}.VY);
            TMP2 = FD_CR4(TMP2);
            
            UXX{k} = TMP1.UX;
            UYY{k} = TMP1.VY;
            VXX{k} = TMP2.UX;
            VYY{k} = TMP2.VY;
            
        end
        
        PG  = cell(length(VEC),1);
        for k = 1:length(VEC)
            PG{k} = struct('x', 0, 'y', 0);
            PG{k}.x = -rho*ACC{k}.x + mu*(UXX{k} + UYY{k}) + rho*BF{k}.x;
            PG{k}.y = -rho*ACC{k}.y + mu*(VXX{k} + VYY{k}) + rho*BF{k}.y;
        end
        
    elseif nargin == 7 && strcmpi(varargin{1}, 'mask');
        
        ACC = PP_TotalAccel(VEC, VGT, dt, 'mask');
        
        UXX = cell(length(VEC),1);
        UYY = cell(length(VEC),1);
        VXX = cell(length(VEC),1);
        VYY = cell(length(VEC),1);
        for k = 1:length(VEC)
            
            TMP1 = struct('C', VEC{k}.C, 'X', VEC{k}.X, 'Y', VEC{k}.Y, ...
                          'U', VGT{k}.UX, 'V', VGT{k}.UY);
            TMP1 = FD_CR4(TMP1, 'mask');
            TMP2 = struct('C', VEC{k}.C, 'X', VEC{k}.X, 'Y', VEC{k}.Y, ...
                          'U', VGT{k}.VX, 'V', VGT{k}.VY);
            TMP2 = FD_CR4(TMP2, 'mask');
            
            UXX{k} = TMP1.UX;
            UYY{k} = TMP1.VY;
            VXX{k} = TMP2.UX;
            VYY{k} = TMP2.VY;
            
        end
        
        PG  = cell(length(VEC),1);
        for k = 1:length(VEC)
            PG{k} = struct('x', 0, 'y', 0);
            PG{k}.x = -rho*ACC{k}.x + mu*(UXX{k} + UYY{k}) + rho*BF{k}.x;
            PG{k}.y = -rho*ACC{k}.y + mu*(VXX{k} + VYY{k}) + rho*BF{k}.y;
        end
        
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
