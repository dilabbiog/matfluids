%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         PREDEFINED FLOW TOOLBOX                         %
%                                                                         %
% doubleGyre                                                              %
% Model Flows                                                             %
% Generate the velocity field of a double gyre                            %
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
% vel = doubleGyre(coord, A, epsi, omega);                                %
% [vel, psi] = doubleGyre(coord, A, epsi, omega);                         %
% [vel, psi, vgt] = doubleGyre(coord, A, epsi, omega);                    %
% [___] = doubleGyre(coord, A, epsi, omega, ctr);                         %
% [___] = doubleGyre(coord, A, epsi, omega, ctr, alpha);                  %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Generate the velocity field for a double gyre [1] with velocity scale   %
% 'A', centre motion amplitude 'epsi' (approximate) and radial frequency  %
% 'omega', where the period is 2*pi/omega. A translation and rotation of  %
% the double gyre can also be specified, the angle being measured         %
% counterclockwise from the positive x axis. The default origin and angle %
% are zero. The stream function and the velocity gradient tensor may also %
% be generated.                                                           %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2019b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
% N/A                                                                     %
%                                                                         %
% Acknowledgments:                                                        %
% N/A                                                                     %
%                                                                         %
% References:                                                             %
% [1] Shadden, S. C., Lekien, F., & Marsden, J. E. (2005). Definition and %
%     properties of Lagrangian coherent structures from finite-time       %
%     Lyapunov exponents in two-dimensional aperiodic flows. Physica D,   %
%     212(3-4), 271-304.                                                  %
%                                                                         %
% ======================================================================= %
% Input Arguments (Required):                                             %
% ----------------------------------------------------------------------- %
% 'coord'        STRUCT ARRAY (1 X 1)                                     %
%              ~ Coordinates. The structure array contains the following  %
%                fields:                                                  %
%                1) coord.x representing the space coordinate x;          %
%                2) coord.y representing the space coordinate y;          %
%                3) coord.t representing the time coordinate t.           %
%                Each field is a column vector (one-dimensional array)    %
%                that strictly increases/decreases monotonically from the %
%                first row to the last.                                   %
% ----------------------------------------------------------------------- %
% 'A'            REAL FINITE SCALAR                                       %
%              ~ Amplitude or velocity scale.                             %
% ----------------------------------------------------------------------- %
% 'epsi'         REAL FINITE SCALAR                                       %
%              ~ Amplitude of the gyre centres (approximate).             %
% ----------------------------------------------------------------------- %
% 'omega'        REAL FINITE SCALAR                                       %
%              ~ Radial frequency of oscillation (period = 2*pi/omga).    %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'ctr'          REAL FINITE TWO-ELEMENT VECTOR                           %
%                Default: [0 0]                                           %
%              ~ Origin location. Translation of the originrelative to    %
%                zero. The origin location is specified as [x y].         %
% ----------------------------------------------------------------------- %
% 'alpha'        REAL FINITE SCALAR                                       %
%                Default: 0                                               %
%              ~ Gyre angle. Rotation of the double gyres by an angle     %
%                measured counterclockwise from the positive x axis.      %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'vel'          STRUCT ARRAY (1 X 1)                                     %
%              ~ Velocity field. The structure array contains the         %
%                following fields:                                        %
%                1) vel.u representing the velocity component along x;    %
%                2) vel.v representing the velocity component along y.    %
%                Each field is a two dimensional array, in ndgrid format  %
%                (i.e., dimensions represent [x y t]).                    %
% ----------------------------------------------------------------------- %
% 'psi'          3-DIMENSIONAL ARRAY                                      %
%              ~ Stream function. The array contains the stream function  %
%                for the double gyre flow. The output array is in ndgrid  %
%                format (i.e., dimensions represent [x y t]).             %
% ----------------------------------------------------------------------- %
% 'vgt'          STRUCT ARRAY (1 X 1)                                     %
%              ~ Velocity gradient tensor. The structure array contains   %
%                the following fields:                                    %
%                1) vgt.ux representing the derivative of u along x;      %
%                2) vgt.uy representing the derivative of u along y;      %
%                3) vgt.vx representing the derivative of v along x;      %
%                4) vgt.vy representing the derivative of v along y.      %
%                Each field is a two dimensional array, in ndgrid format  %
%                (i.e., dimensions represent [x y t]).                    %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Generate a double gyre on the domain (x,y) = ([0,2],[0,1]) with a       %
% constant grid spacing of 0.01 over the time interval [0,20] with time   %
% step size 0.1. Use A = 0.1, epsi = 0.25, and omega = 2*pi/10.           %
%                                                                         %
% >> coord.x = linspace(0, 2, 201).';                                     %
% >> coord.y = linspace(0, 1, 101).';                                     %
% >> coord.t = linspace(0, 20, 21).';                                     %
% >> A       = 0.1;                                                       %
% >> epsi    = 0.25;                                                      %
% >> omega   = 2*pi/10;                                                   %
% >> vel     = doubleGyre(coord, A, epsi, omega);                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [vel, varargout] = doubleGyre(coord, A, epsi, omega, varargin)


%% PARSE INPUTS

% Input defaults.
default.ctr   = [0 0];
default.alpha = 0;

% Input checks.
check.coord = @(x) examineCoord(x, ["x" "y" "t"], 'all');
check.A     = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'finite', 'real', 'scalar'});
check.epsi  = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'finite', 'real', 'scalar'});
check.omega = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'finite', 'real', 'scalar'});
check.ctr   = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'finite', 'real', 'vector', 'numel', 2});
check.alpha = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'finite', 'real', 'scalar'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'coord' ,                 check.coord );
addRequired ( hParser, 'A'     ,                 check.A     );
addRequired ( hParser, 'epsi'  ,                 check.epsi  );
addRequired ( hParser, 'omega' ,                 check.omega );
addOptional ( hParser, 'ctr'   , default.ctr   , check.ctr   );
addOptional ( hParser, 'alpha' , default.alpha , check.alpha );
parse(hParser, coord, A, epsi, omega, varargin{:});
clear check default;

% Additional verifications.
narginchk(4,6);
nargoutchk(0,3);


%% STREAM FUNCTION

% Transform the coordinates.
ctr     = hParser.Results.ctr;
alpha   = hParser.Results.alpha;
[X0,Y0] = ndgrid(coord.x, coord.y);
X       =  (X0 - ctr(1))*cos(alpha) + (Y0 - ctr(2))*sin(alpha);
Y       = -(X0 - ctr(1))*sin(alpha) + (Y0 - ctr(2))*cos(alpha);
clear ctr X0 Y0;

% Define intermediate variables.
A     = hParser.Results.A;
epsi  = hParser.Results.epsi;
omega = hParser.Results.omega;
t0    = permute(coord.t, [3 2 1]);
at    = epsi*sin(omega*t0);
bt    = 1 - 2*at;
fxt   = at.*X.^2 + bt.*X;
dfxt  = 2*at.*X + bt;
d2fxt = 2*at;
clear at bt t0;

% Define the stream function.
psi = -A*sin(pi*fxt).*sin(pi*Y);
if nargout > 1, varargout{1} = psi; end
clear psi;


%% VELOCITY FIELD

% Initialize the velocity field.
vel = struct('u' , 0, 'v' , 0);

% Compute the velocity field.
vel.u = -A*pi*cos(alpha)*sin(pi*fxt).*cos(pi*Y);
vel.v =  A*pi*cos(alpha)*dfxt.*cos(pi*fxt).*sin(pi*Y);


%% VELOCITY GRADIENT TENSOR

if nargout > 2
    % Initialize the velocity gradient tensor.
    varargout{2} = struct('ux' , 0, 'uy' , 0, 'vx' , 0, 'vy' , 0);

    % Compute the velocity gradient tensor.
    varargout{2}.ux = -A*pi^2*cos(alpha)*dfxt.*cos(pi*fxt).*cos(pi*Y);
    varargout{2}.uy =  A*pi^2*cos(alpha)*sin(pi*fxt).*sin(pi*Y);
    varargout{2}.vx =  A*pi*cos(alpha)*sin(pi*Y)                        ...
                   .*( d2fxt.*cos(pi*fxt)                               ...
                    -  pi*dfxt.^2.*sin(pi*fxt));
    varargout{2}.vy =  A*pi^2*cos(alpha)*dfxt.*cos(pi*fxt).*cos(pi*Y);
end
clear dxt dfxt d2fxt X Y;


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
% 2022/06/15 -- (GDL) Fixed input check for coord.                        %
% 2022/06/10 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
