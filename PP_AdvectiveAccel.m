%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                   ADVECTIVE ACCELERATION VECTOR FIELD                   %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Giuseppe Di Labbio                                                      %
% Department of Mechanical, Industrial & Aerospace Engineering            %
% Concordia University Montréal, Canada                                   %
%                                                                         %
% Last Update: October 2nd, 2018 by Giuseppe Di Labbio                    %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% SYNTAX                                                                  %
%                                                                         %
% ACC = PP_AdvectiveAccel(VEC, VGT);                                      %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the advective acceleration vector field of a     %
% velocity field given its velocity gradient tensor. The function detects %
% if the calculation is two or three dimensional.                         %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'ACC'        - STRUCT                                                   %
%              - Advective acceleration vector field in the form of a     %
%                struct. The struct contains the vector components of the %
%                acceleration field ('x' and 'y' in two dimensions). In   %
%                three dimensions, the struct ought to contain the        %
%                components 'x', 'y' and 'z'.                             %
% ----------------------------------------------------------------------- %
% 'VEC'        - STRUCT                                                   %
%              - Velocity vector field in the form of a struct. The       %
%                struct contains the mask ('C'), Cartesian coordinates    %
%                ('X' and 'Y') and velocity vector components ('U' and    %
%                'V'). In three dimensions, the struct ought to contain   %
%                the Cartesian coordinates 'X', 'Y', and 'Z' and the      %
%                velocity components 'U', 'V', and 'W'.                   %
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
% Calculate the advective acceleration vector field of a time-dependent   %
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
% >> ACC  = cell(length(t),1);                                            %
% >> Amag = cell(length(t),1);                                            %
% >> for k = 1:length(t)                                                  %
% ACC{k} = PP_AdvectiveAccel(VEC{k}, VGT{k});                             %
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
% N/A                                                                     %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_AdvectiveAccel
function [ACC] = PP_AdvectiveAccel(VEC, VGT)

if length(fieldnames(VEC)) == 5
    
    ACC = struct('x', 0, 'y', 0);
    
    ACC.x = VEC.U.*VGT.UX + VEC.V.*VGT.UY;
    ACC.y = VEC.U.*VGT.VX + VEC.V.*VGT.VY;
    
elseif length(fieldnames(VEC)) == 7
    
    ACC = struct('x', 0, 'y', 0, 'z', 0);

    ACC.x = VEC.U.*VGT.UX + VEC.V.*VGT.UY + VEC.W.*VGT.UZ;
    ACC.y = VEC.U.*VGT.VX + VEC.V.*VGT.VY + VEC.W.*VGT.VZ;
    ACC.z = VEC.U.*VGT.WX + VEC.V.*VGT.WY + VEC.W.*VGT.WZ;
    
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
