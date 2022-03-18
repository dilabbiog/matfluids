%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            ADVECTION TOOLBOX                            %
%                                                                         %
% advectionODE                                                            %
% Particle Advection                                                      %
% Basic ordinary differential equation used for particle advection        %
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
% dxdt = advectionODE(t, x, coord, vel, method, extrapval);               %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Basic ordinary differential equation to use for particle advection. The %
% function interpolates the velocity field onto the query points.         %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2021b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
% advectionODE                                                            %
%                                                                         %
% Acknowledgments:                                                        %
% N/A                                                                     %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ======================================================================= %
% Input Arguments (Required):                                             %
% ----------------------------------------------------------------------- %
% 't'            REAL SCALAR                                              %
%              ~ Current solution time.                                   %
% ----------------------------------------------------------------------- %
% 'x'            REAL ARRAY                                               %
%              ~ Particle positions at the current solution time.         %
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
% 'method'       CASE-INSENSITIVE CHARACTER ARRAY                         %
%              ~ Specify the method of interpolation for intermediate     %
%                grid points. See the MATLAB documentation of interpn for %
%                more information.                                        %
% ----------------------------------------------------------------------- %
% 'extrapval'    SCALAR                                                   %
%              ~ Scalar value assigned to query points falling outside    %
%                the domain. Unlike MATLAB's default for interpn, where   %
%                extrapolation is performed for the 'makima' and 'spline' %
%                methods and NaN is assigned for others, here only scalar %
%                assignment is permitted, regardless of the interpolation %
%                method.                                                  %
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
% 'dxdt'         REAL ARRAY                                               %
%              ~ Velocity of the particles at the current solution time.  %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dxdt] = advectionODE(t, x, coord, vel, method, extrapval)


%% PARSE INPUTS

% Input defaults.
% N/A

% Input checks.
% N/A

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 't'         );
addRequired ( hParser, 'x'         );
addRequired ( hParser, 'coord'     );
addRequired ( hParser, 'vel'       );
addRequired ( hParser, 'method'    );
addRequired ( hParser, 'extrapval' );
parse(hParser, t, x, coord, vel, method, extrapval);

% Additional verifications.
narginchk(6,6);
nargoutchk(0,1);


%% ADVECTION ODE

% Get the fields in the coordinate and velocity field structures.
fieldsCoord    = fieldnames(coord);
numFieldsCoord = length(fieldsCoord);
fieldsVel      = fieldnames(vel);

% Reshape the query points.
n = numel(x)/(numFieldsCoord - 1);
x = reshape(x, [], n);

% Create a cell array to hold the coordinates.
coordCell = cell(numFieldsCoord,1);
for k = 1:numFieldsCoord
    coordCell{k} = coord.(fieldsCoord{k});
end

% Create a cell array to hold the query points.
xCell = cell(numFieldsCoord,1);
for k = 1:numFieldsCoord-1
    xCell{k} = x(k,:);
end
xCell{numFieldsCoord} = t*ones(1,n);

% Interpolate the velocity field at the query points.
dxdt = zeros(numFieldsCoord-1, n);
for k = 1:numFieldsCoord-1
    dxdt(k,:) = interpn(coordCell{:}, vel.(fieldsVel{k}), xCell{:},     ...
                        method, extrapval);
end
dxdt = dxdt(:);


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
% 2022/03/18 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% * Limit the time interval used for the interpolation (to be tested).    %
%   [~,idxt] = min(abs(coord.t - t));                                     %
%   idx = repmat({':'}, 1, numFieldsCoord-1);                             %
%   nt = length(coord.t);                                                 %
%   idxt = idxt + [-(idxt == nt); (idxt < nt)];                           %
%   coordCell{numFieldsCoord} = coord.t(idxt);                            %
%   vel.(fieldsVel{k}) --> vel.(fieldsVel{k})(idx{:},idxt);               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
