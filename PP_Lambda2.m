%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                            LAMBDA-2 CRITERION                           %
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
% L2 = PP_Lambda2(VGT);                                                   %
% L2 = PP_Lambda2(VGT, opts);                                             %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the Lambda2-criterion for 2D or 3D data as       %
% defined by Jeong and Hussain [1]. The user can optionally choose to     %
% keep only negative values of Lambda2 and to normalize the values of     %
% Lambda2 with respect to its largest negative value (in absolute).       %
% Recall that negative values of Lambda2 signify the presence of a        %
% vortex.                                                                 %
%                                                                         %
% References:                                                             %
% [1] Jeong, J. & Hussain, F. (1995). On  the identification of a vortex. %
%     Journal of Fluid Mechanics, 285, 69-94.                             %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'L2'         - 2D DOUBLE ARRAY                                          %
%              - Lambda2 criterion.                                       %
% ----------------------------------------------------------------------- %
% 'opts'       - SPECIFIC STRING                                          %
%              - Option to normalize the values of lambda-2 by its        %
%                largest negative value (in absolute) using 'norm' and to %
%                retain only negative values of lambda-2 using 'neg'.     %
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
% Calculate the Lambda2 criterion of a time-dependent double gyre on the  %
% domain (x,y) = [0,2]x[0,1] with a constant grid spacing of 0.01 over    %
% the time interval [0,20] with time-step size 0.1. Use A = 0.1, epsilon  %
% = 0.25, and omega = 2*pi/10. Keep only negative values of normalized    %
% Lambda2.                                                                %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> L2 = cell(length(t),1);                                              %
% >> for k = 1:length(t)                                                  %
% L2{k} = PP_Lambda2(VGT{k}, 'neg', 'norm');                              %
% end                                                                     %
% >> for k = 1:length(t)                                                  %
% contourf(VEC{k}.X, VEC{k}.Y, L2{k}, 'EdgeColor', 'None');               %
% colormap gray;                                                          %
% hold on;                                                                %
% quiver(VEC{k}.X(1:4:end,1:4:end), VEC{k}.Y(1:4:end,1:4:end), ...        %
%        VEC{k}.U(1:4:end,1:4:end), VEC{k}.V(1:4:end,1:4:end), 'b');      %
% hold off;                                                               %
% axis([0 2 0 1]);                                                        %
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

%% PP_Lambda2
function [L2] = PP_Lambda2(VGT, varargin)

if length(fieldnames(VGT)) == 4
    
    [ny, nx] = size(VGT.UX);
    L2 = zeros(ny, nx);
    
    M = zeros(2);
    for i = 1:ny
        for j = 1:nx
            
            M(1,1) = 2*VGT.UX(i,j)^2 + 2*VGT.UY(i,j)*VGT.VX(i,j);
            M(1,2) = VGT.UX(i,j)*VGT.UY(i,j) + VGT.UY(i,j)*VGT.VY(i,j) ...
                   + VGT.UX(i,j)*VGT.VX(i,j) + VGT.VX(i,j)*VGT.VY(i,j);
            M(2,1) = M(1,2);
            M(2,2) = 2*VGT.VY(i,j)^2 + 2*VGT.UY(i,j)*VGT.VX(i,j);
            
            L2(i,j) = min(eig(0.5*M));
            
        end
    end
    
    switch nargin
        case 2
            if strcmpi(varargin{1},'neg')
                L2(L2 > 0) = 0;
            end
            
            if strcmpi(varargin{1},'norm')
                L2 = L2/max(max(-L2));
            end
        case 3
            if strcmpi(varargin{1},'neg') || strcmpi(varargin{2},'neg')
                L2(L2 > 0) = 0;
            end
            
            if strcmpi(varargin{1},'norm') || strcmpi(varargin{2},'norm')
                L2 = L2/max(max(-L2));
            end
    end
    
elseif length(fieldnames(VGT)) == 9
    
    [ny, nx, nz] = size(VGT.UX);
    L2 = zeros(ny, nx, nz);
    
    M = zeros(3);
    for i = 1:ny
        for j = 1:nx
            for k = 1:nz
            
                M(1,1) = 2*VGT.UX(i,j)^2 + 2*VGT.UY(i,j)*VGT.VX(i,j) ...
                       + 2*VGT.UZ(i,j)*VGT.WX(i,j);
                M(1,2) = VGT.UX(i,j)*VGT.UY(i,j) ...
                       + VGT.UY(i,j)*VGT.VY(i,j) ...
                       + VGT.UZ(i,j)*VGT.WY(i,j) ...
                       + VGT.UX(i,j)*VGT.VX(i,j) ...
                       + VGT.VX(i,j)*VGT.VY(i,j) ...
                       + VGT.WX(i,j)*VGT.VZ(i,j);
                M(1,3) = VGT.UX(i,j)*VGT.UZ(i,j) ...
                       + VGT.UY(i,j)*VGT.VZ(i,j) ...
                       + VGT.UZ(i,j)*VGT.WZ(i,j) ...
                       + VGT.UX(i,j)*VGT.WX(i,j) ...
                       + VGT.VX(i,j)*VGT.WY(i,j) ...
                       + VGT.WX(i,j)*VGT.WZ(i,j);
                M(2,1) = M(1,2);
                M(2,2) = 2*VGT.VY(i,j)^2 + 2*VGT.UY(i,j)*VGT.VX(i,j) ...
                       + 2*VGT.VZ(i,j)*VGT.WY(i,j);
                M(2,3) = VGT.VX(i,j)*VGT.UZ(i,j) ...
                       + VGT.VY(i,j)*VGT.VZ(i,j) ...
                       + VGT.VZ(i,j)*VGT.WZ(i,j) ...
                       + VGT.UY(i,j)*VGT.WX(i,j) ...
                       + VGT.VY(i,j)*VGT.WY(i,j) ...
                       + VGT.WY(i,j)*VGT.WZ(i,j);
                M(3,1) = M(1,3);
                M(3,2) = M(2,3);
                M(3,3) = 2*VGT.WZ(i,j)^2 + 2*VGT.UZ(i,j)*VGT.WX(i,j) ...
                       + 2*VGT.VZ(i,j)*VGT.WY(i,j);
                
                E = sort(eig(0.5*M));
                L2(i,j,k) = E(2);
            
            end
        end
    end
    
    switch nargin
        case 2
            if strcmpi(varargin{1},'neg')
                L2(L2 > 0) = 0;
            end
            
            if strcmpi(varargin{1},'norm')
                L2 = L2/max(max(max(-L2)));
            end
        case 3
            if strcmpi(varargin{1},'neg') || strcmpi(varargin{2},'neg')
                L2(L2 > 0) = 0;
            end
            
            if strcmpi(varargin{1},'norm') || strcmpi(varargin{2},'norm')
                L2 = L2/max(max(max(-L2)));
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
