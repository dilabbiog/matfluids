%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             VORTEX TOOLBOX                              %
%                                                                         %
% gamma2                                                                  %
% Vortex Identification                                                   %
% Compute the Gamma_2 criterion to identify vortex cores                  %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Département de génie mécanique                                          %
% École de technologie supérieure (ÉTS)                                   %
% Montréal, Québec                                                        %
% Canada                                                                  %
%                                                                         %
% Contributors: Giuseppe Di Labbio                                        %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Copyright (C) 2022 Giuseppe Di Labbio                                   %
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
% G2 = gamma2(coord, vel, winSize);                                       %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Compute the Gamma_2 criterion using a specified window size to identify %
% vortex cores (cf. [1-2]) in a two-dimensional, time-dependent velocity  %
% field.                                                                  %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2021b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
% N/A                                                                     %
%                                                                         %
% Acknowledgments:                                                        %
% [1] I am not the first to realize that MATLAB's convolution functions   %
%     can be used to compute the Gamma_1 and Gamma_2 criteria. See for    %
%     example Fernando Zigunov's implementation of the Gamma_1 criterion  %
%     using the conv2 function on MATLAB's File Exchange. Kudos!          %
%                                                                         %
% References:                                                             %
% [1] Michard, M., Graftieaux, L., Lollini, L., & Grosjean, N. (1997).    %
%     Identification of vortical structures by a non local criterion:     %
%     Application to PIV measurements and DNS-LES results of turbulent    %
%     rotating flows. In Proceedings of the 11th Symposium on Turbulent   %
%     Shear Flows -- TSF11 (Vol. 3, pp. 28-25--28-30). Grenoble, FR:      %
%     Institut National Polytechnique de Grenoble.                        %
% [2] Graftieaux, L., Michard, M., & Grosjean, N. (2001). Combining PIV,  %
%     POD and vortex identification algorithms for the study of unsteady  %
%     turbulent swirling flows. Measurement Science & Technology, 12(9),  %
%     1422--1429.                                                         %
%                                                                         %
% ======================================================================= %
% Input Arguments (Required):                                             %
% ----------------------------------------------------------------------- %
% 'coord'        STRUCT ARRAY (1 X 1)                                     %
%              ~ Coordinates. The structure array may contain the         %
%                following fields:                                        %
%                1) coord.x representing the space coordinate x;          %
%                2) coord.y representing the space coordinate y;          %
%                3) coord.z representing the space coordinate z;          %
%                4) coord.t representing the time coordinate t.           %
%                Each field is a column vector (one-dimensional array)    %
%                that strictly increases/decreases monotonically from the %
%                first row to the last. Not all fields need be present in %
%                the coordinate structure, however the order must follow  %
%                [x y z t]. For example, the coordinate structure may     %
%                contain fields ordered as [x z t] but not as [t x z].    %
%                Both empty and scalar arrays are permitted.              %
% ----------------------------------------------------------------------- %
% 'vel'          STRUCT ARRAY (1 X 1)                                     %
%              ~ Velocity field. The structure array may contain the      %
%                following fields:                                        %
%                1) vel.u representing the x-component velocity field;    %
%                2) vel.v representing the y-component velocity field;    %
%                3) vel.w representing the z-component velocity field.    %
%                Each field is a two- to four-dimensional array, with the %
%                dimensions corresponding to the order [x y z t]. Not all %
%                fields need be present in the velocity field structure,  %
%                however the order must follow [u v w]. For example, the  %
%                velocity field structure may contain fields ordered as   %
%                [u w] but not as [w u]. Both empty and scalar arrays are %
%                permitted.                                               %
% ----------------------------------------------------------------------- %
% 'winSize'      POSITIVE INTEGER SCALAR                                  %
%              ~ Kernel window size. Size of the sliding sample window to %
%                use for computing the Gamma_2 criterion. Ideally an odd  %
%                integer, but no restriction is placed on parity.         %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'G2'           REAL ARRAY                                               %
%              ~ Gamma_2 criterion. The output will have the same size as %
%                the input velocity field. A value of closer to +1        %
%                indicates more coherent counterclockwise rotation with a %
%                value closer to -1 indicates more coherent clockwise     %
%                rotation.                                                %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Compute the Gamma_2 criterion of a double gyre. Generate the double     %
% gyre on the domain (x,y) = ([0,2],[0,1]) using a constant grid spacing  %
% of 0.005 over the time interval [0,20] with time step size 0.1. Use     %
% A = 0.1, epsi = 0.25, and omega = 2*pi/10.                              %
%                                                                         %
% >> coord.x = linspace(0, 2, 401).';                                     %
% >> coord.y = linspace(0, 1, 201).';                                     %
% >> coord.t = linspace(0, 20, 201).';                                    %
% >> A       = 0.1;                                                       %
% >> epsi    = 0.25;                                                      %
% >> omega   = 2*pi/10;                                                   %
% >> vel     = doubleGyre(coord, A, epsi, omega);                         %
% >> G2      = gamma2(coord, vel, 11);                                    %
% >> pcolor(coord.x, coord.y, G2(:,:,26).');                              %
% >> shading interp;                                                      %
% >> axis equal tight;                                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [G2] = gamma2(coord, vel, winSize)


%% PARSE INPUTS

% Input defaults.
% N/A

% Input checks.
check.winSize = @(x) validateattributes(x,                              ...
                     {'logical', 'numeric'},                            ...
                     {'finite', 'positive', 'real'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'coord'                  );
addRequired ( hParser, 'vel'                    );
addRequired ( hParser, 'winSize', check.winSize );
parse(hParser, coord, vel, winSize);
clear check default;

% Additional verifications.
narginchk(3,3);
nargoutchk(0,1);


%% COMPUTE GAMMA_2

% Determine the spatial discretizations.
d.x = diff(coord.x(1:2));
d.y = diff(coord.y(1:2));

% Define the kernels.
% Recall that the kernels will be mirrored by convn.
kernAvg = ones(winSize);
kernPMx = d.x*repmat((winSize-1)/2:-1:-(winSize-1)/2, [winSize 1]).';
kernPMy = d.y*repmat((winSize-1)/2:-1:-(winSize-1)/2, [winSize 1]);

% Normalize the cross product kernel.
magPM   = sqrt(kernPMx.^2 + kernPMy.^2);
kernPMx = kernPMx./magPM;
kernPMy = kernPMy./magPM;
kernPMx(isnan(kernPMx)) = 0;
kernPMy(isnan(kernPMy)) = 0;

% Compute Ump = Um - Up and Vmp = Vm - Vp.
Ump = vel.u - convn(vel.u, kernAvg, 'same')/numel(kernAvg);
Vmp = vel.v - convn(vel.v, kernAvg, 'same')/numel(kernAvg);

% Normalize Ump and Vmp.
magVel = sqrt(Ump.^2 + Vmp.^2);
Ump    = Ump./magVel;
Vmp    = Vmp./magVel;
Ump(isnan(Ump)) = 0;
Vmp(isnan(Vmp)) = 0;

% Determine the area of the window.
winArea = winSize.^2;

% Compute Gamma_2.
G2 = (convn(Vmp, kernPMx, 'same') - convn(Ump, kernPMy, 'same'))/winArea;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  NOTES                                  %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Line(s) N/A                                                             %
% * N/A                                                                   %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% SUPPRESSED MESSAGES                                                     %
%                                                                         %
% Line(s) N/A                                                             %
% Message(s)                                                              %
% * N/A                                                                   %
% Reason(s)                                                               %
% * N/A                                                                   %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% CHANGE LOG                                                              %
%                                                                         %
% 2022/06/17 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
