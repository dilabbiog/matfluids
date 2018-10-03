%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                     CAUCHY-GREEN DEFORMATION TENSOR                     %
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
% CG = PP_CauchyGreenTensor(A, ds, scheme, dir);                          %
% CG = PP_CauchyGreenTensor(A, ds, scheme, dir, mask);                    %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the left or right Cauchy-Green deformation       %
% gradient tensor, also known as the Finger and Cauchy tensors, of an     %
% advected particle field. The user may specify the finite difference     %
% scheme to be used and whether or not to use a pre-defined mask.         %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'A'          - STRUCT                                                   %
%              - Struct containing fields X and Y holding the positions   %
%                of the advected particles.                               %
% ----------------------------------------------------------------------- %
% 'dir'        - SPECIFIC STRING                                          %
%              - String specifying whether to compute the 'left' or       %
%                'right' Cauchy-Green deformation tensor.                 %
% ----------------------------------------------------------------------- %
% 'ds'         - TWO-ELEMENT REAL VECTOR                                  %
%              - Spatial grid spacing in the x and y directions. The user %
%                must specify ds in the form ds = [dx dy] (or             %
%                equivalently ds = [dx; dy]).                             %
% ----------------------------------------------------------------------- %
% 'mask'      - 2D DOUBLE ARRAY                                           %
%             - Mask information for the vector field. Values of 1        %
%               correspond to active nodes while values of 0 correspond   %
%               to blanked nodes.                                         %
% ----------------------------------------------------------------------- %
% 'scheme'    - SPECIFIC STRING                                           %
%             - String denoting which scheme to use. The user can select  %
%               from the explicit second order central scheme ('EXO2'),   %
%               or Lele's second, third or fourth order compact schemes   %
%               ('COMP2', 'COMP3', 'COMP4').                              %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Advect a grid of points placed at each node of a time-dependent double  %
% gyre over the domain (x,y) = [0,2]x[0,1] with a constant grid spacing   %
% of 0.01 over the time interval [0,20] with time-step size 0.1. Use A =  %
% 0.1, epsilon = 0.25, and omega = 2*pi/10. Use the time-step dt = 0.1    %
% for the advection (i.e. use the 'singlestep' option). Compute the right %
% Cauchy-Green deformation tensor of the advected particle field at the   %
% final advection time using an explicit second-order centered finite     %
% difference scheme.                                                      %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> ADV = PP_AdvectGrid(VEC{1}, VEC, t(2), 'singlestep');                %
% >> CG = PP_CauchyGreenTensor(ADV{21}, [x(2) y(2)], 'EXO2', 'right');    %
% >> contourf(VEC{21}.X, VEC{21}.Y, CG.yx, 'EdgeColor', 'None');          %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% FD_COMP2                                                                %
% FD_COMP3                                                                %
% FD_COMP4                                                                %
% FD_EXO2                                                                 %
% FD_SCH                                                                  %
% PP_DeformGradTensor                                                     %
% TDMA                                                                    %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_CauchyGreenTensor
function [CG] = PP_CauchyGreenTensor(A, ds, scheme, dir, varargin)

% Compute the deformation gradient tensor.
if nargin == 5
    F = PP_DeformGradTensor(A, ds, scheme, varargin{1});
else
    F = PP_DeformGradTensor(A, ds, scheme);
end

[ny,nx] = size(F.xx);
CG = struct('xx', zeros(ny,nx), 'xy', zeros(ny,nx), ...
            'yx', zeros(ny,nx), 'yy', zeros(ny,nx));
if strcmpi(dir, 'left')
    for i = 1:ny
        for j = 1:nx
            
            M = [F.xx(i,j) F.xy(i,j); F.yx(i,j) F.yy(i,j)];
            M = M*M';
            CG.xx(i,j) = M(1,1);
            CG.xy(i,j) = M(1,2);
            CG.yx(i,j) = M(2,1);
            CG.yy(i,j) = M(2,2);
            
        end
    end
elseif strcmpi(dir, 'right')
    for i = 1:ny
        for j = 1:nx
            
            M = [F.xx(i,j) F.xy(i,j); F.yx(i,j) F.yy(i,j)];
            M = M'*M;
            CG.xx(i,j) = M(1,1);
            CG.xy(i,j) = M(1,2);
            CG.yx(i,j) = M(2,1);
            CG.yy(i,j) = M(2,2);
            
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
