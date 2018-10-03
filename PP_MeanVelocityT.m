%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                          MEAN TEMPORAL VELOCITY                         %
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
% Vmean = PP_MeanVelocityT(VEC);                                          %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the mean temporal velocity of a velocity field.  %
% The function detects if the calculation is two or three dimensional.    %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
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
% 'Vmean'      - 2D DOUBLE ARRAY                                          %
%              - Mean temporal velocity calculated in two or three        %
%                dimensions.                                              %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Calculate the mean temporal velocity of a time-dependent double gyre on %
% the domain (x,y) = [0,2]x[0,1], with a constant grid spacing of 0.01,   %
% over the time interval [0,20] with time-step size 0.1. Use A = 0.1,     %
% epsilon = 0.25 and omega = 2*pi/10.                                     %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> Vmean = PP_MeanVelocityT(VEC);                                       %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% field2mat3D                                                             %
% field2mat4D                                                             %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_MeanVelocityT
function [Vmean] = PP_MeanVelocityT(VEC)

if length(fieldnames(VEC{1})) == 5
    %% Two dimensional case.
    
    Vmean = VEC{1};
    Vmean.U = mean(field2mat3D(VEC, 'U'), 3);
    Vmean.V = mean(field2mat3D(VEC, 'V'), 3);
    
elseif length(fieldnames(VEC{1})) == 7
    %% Three dimensional case.
    
    Vmean = VEC{1};
    Vmean.U = mean(field2mat4D(VEC, 'U'), 4);
    Vmean.V = mean(field2mat4D(VEC, 'V'), 4);
    Vmean.W = mean(field2mat4D(VEC, 'W'), 4);
    
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
