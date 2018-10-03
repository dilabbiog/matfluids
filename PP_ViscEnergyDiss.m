%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                     VISCOUS ENERGY DISSIPATION RATE                     %
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
% VED = ViscEnergyDiss(VGT, mu);                                          %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the viscous eneregy dissipation rate of a        %
% velocity field given its velocity gradient tensor and dynamic           %
% viscosity. The function detects if the calculation is two or three      %
% dimensional.                                                            %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'mu'         - REAL SCALAR                                              %
%              - Dynamic viscosity of the fluid.                          %
% ----------------------------------------------------------------------- %
% 'VED'        - 2D DOUBLE ARRAY                                          %
%              - Viscous energy dissipation rate.                         %
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
% Calculate the viscous energy dissipation rate of a time-dependent       %
% double gyre on the domain (x,y) = [0,2]x[0,1] with a constant grid      %
% spacing of 0.01 over the time interval [0,20] with time-step size 0.1.  %
% Use A = 0.1, epsilon = 0.25, omega = 2*pi/10, and mu = 0.001.           %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> mu   = 0.001;                                                        %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> VED = cell(length(t),1);                                             %
% >> for k = 1:length(t)                                                  %
% VED{k} = PP_ViscEnergyDiss(VGT{k}, mu);                                 %
% end                                                                     %
% >> for k = 1:length(t)                                                  %
% contourf(VEC{k}.X, VEC{k}.Y, VED{k}, 'EdgeColor', 'None');              %
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

%% PP_ViscEnergyDiss
function [VED] = PP_ViscEnergyDiss(VGT, mu)

if length(fieldnames(VGT)) == 4
    
    VED = (mu/2)*((2*VGT.UX).^2 + 2*(VGT.UY + VGT.VX).^2 + (2*VGT.VY).^2);
    
elseif length(fieldnames(VGT)) == 9
    
    VED = (mu/2)*(2*(VGT.UY + VGT.VX).^2 + (2*VGT.UX).^2 ...
        +         2*(VGT.UZ + VGT.WX).^2 + (2*VGT.VY).^2 ...
        +         2*(VGT.VZ + VGT.WY).^2 + (2*VGT.WZ).^2);
   
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
