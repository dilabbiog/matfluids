%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                         NEWTONIAN STRESS TENSOR                         %
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
% NST = PP_NewtStressTensor(VGT, mu);                                     %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the Newtonian stress tensor of a velocity field  %
% given its velocity gradient tensor and dynamic viscosity. The function  %
% detects if the calculation is two or three dimensional.                 %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'mu'         - REAL SCALAR                                              %
%              - Dynamic viscosity.                                       %
% ----------------------------------------------------------------------- %
% 'NST'        - STRUCT                                                   %
%              - Newtonian stress tensor containing the four components   %
%                'xx', 'xy', 'yx' and 'yy' in 2D or the nine components   %
%                'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy' and 'zz'  %
%                in 3D.                                                   %
% ----------------------------------------------------------------------- %
% 'VGT'        - STRUCT                                                   %
%              - Velocity gradient tensor containing the four components  %
%                'UX', 'UY', 'VX' and 'VY' in 2D or the nine components   %
%                'UX', 'UY', 'UZ', 'VX', 'VY', 'VZ', 'WX', 'WY' and 'WZ'  %
%                in 3D.                                                   %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Calculate the Newtonian stress tensor of a time-dependent double gyre   %
% on the domain (x,y) = [0,2]x[0,1] with a constant grid spacing of 0.01  %
% over the time interval [0,20] with time-step size 0.1. Use A = 0.1,     %
% epsilon = 0.25, omega = 2*pi/10, and mu = 0.001.                        %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> mu   = 0.001;                                                        %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> NST = cell(length(t),1);                                             %
% >> for k = 1:length(t)                                                  %
% NST{k} = PP_NewtStressTensor(VGT{k}, mu);                               %
% end                                                                     %
% >> for k = 1:length(t)                                                  %
% contourf(VEC{k}.X, VEC{k}.Y, NST{k}.yx, 'EdgeColor', 'None');           %
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
% N/A                                                                     %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_NewtStressTensor
function [NST] = PP_NewtStressTensor(VGT, mu)

if length(fieldnames(VGT)) == 4
    
    NST = struct('xx', 0, 'xy', 0, 'yx', 0, 'yy', 0);
    
    NST.xx = mu*(VGT.UX + VGT.UX);
    NST.xy = mu*(VGT.UY + VGT.VX);
    NST.yx = NST.xy;
    NST.yy = mu*(VGT.VY + VGT.VY);
    
elseif length(fieldnames(VGT)) == 9
    
    NST = struct('xx', 0, 'xy', 0, 'xz', 0, 'yx', 0, 'yy', 0, 'yz', 0, ...
                 'zx', 0, 'zy', 0, 'zz', 0);
    
    NST.xx = mu*(VGT.UX + VGT.UX);
    NST.xy = mu*(VGT.UY + VGT.VX);
    NST.xz = mu*(VGT.UZ + VGT.WX);
    NST.yx = NST.xy;
    NST.yy = mu*(VGT.VY + VGT.VY);
    NST.yz = mu*(VGT.VZ + VGT.WY);
    NST.zx = NST.xz;
    NST.zy = NST.yz;
    NST.zz = mu*(VGT.WZ + VGT.WZ);
    
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
