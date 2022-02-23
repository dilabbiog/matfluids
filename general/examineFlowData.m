%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% examineFlowData                                                         %
% MATfluids Structure Interrogation                                       %
% Examine contents of both the coordinate and velocity field structures   %
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
% examineFlowData(coord, vel);                                            %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Examine whether the coordinate and velocity field structures are        %
% properly formatted for MATfluids. The function has no specific output,  %
% it will simply throw an error if either the coordinate or velocity      %
% field structures is not correctly formatted, or if their sizes are      %
% incompatible.                                                           %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2019b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
% examineCoord                                                            %
% examineVel                                                              %
% isnegreal                                                               %
% isposreal                                                               %
% isrealnum                                                               %
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
%                1) vel.u representing the velocity component along x;    %
%                2) vel.v representing the velocity component along y;    %
%                3) vel.w representing the velocity component along z.    %
%                Each field is a one to four dimensional array, in ndgrid %
%                format (i.e., dimensions represent [x y z t]). All       %
%                constituent arrays must have the same size. Not all      %
%                fields need be present in the velocity field structure,  %
%                however the order must follow [u v w]. For example, the  %
%                velocity field structure may contain fields ordered as   %
%                [u w] but not as [w u]. Both empty and scalar arrays are %
%                permitted.                                               %
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
% N/A                                                                     %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Determine whether the following coordinate and velocity field           %
% structures are correct for MATfluids.                                   %
%                                                                         %
% >> clear coord vel;                                                     %
% >> coord.x = (-10:0.1:10).';                                            %
% >> coord.y = (0:0.1:5).';                                               %
% >> coord.t = (0:0.1:2).';                                               %
% >> vel.u   = 2;                                                         %
% >> vel.v   = rand(length(coord.x), length(coord.y), length(coord.t));   %
% >> examineFlowData(coord, vel);                                         %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Determine whether the following coordinate and velocity field           %
% structures are correct for MATfluids.                                   %
%                                                                         %
% >> clear coord vel;                                                     %
% >> coord.x = (-2:0.05:2).';                                             %
% >> coord.z = (0:0.1:1).';                                               %
% >> coord.t = (0:0.01:1).';                                              %
% >> vel.u   = rand(length(coord.x), length(coord.z), length(coord.t));   %
% >> vel.v   = [];                                                        %
% >> vel.w   = rand(length(coord.x), length(coord.z), length(coord.t));   %
% >> examineFlowData(coord, vel);                                         %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 3                                                               %
%                                                                         %
% Determine whether the following coordinate and velocity field           %
% structures are correct for MATfluids.                                   %
%                                                                         %
% >> clear coord vel;                                                     %
% >> coord.x = (-2:0.05:2).';                                             %
% >> coord.y = 1;                                                         %
% >> coord.z = (0:0.1:1).';                                               %
% >> coord.t = (0:0.01:1).';                                              %
% >> vel.u   = rand(length(coord.x), length(coord.z), length(coord.t));   %
% >> vel.w   = rand(length(coord.x), length(coord.z), length(coord.t));   %
% >> examineFlowData(coord, vel);                                         %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 4                                                               %
%                                                                         %
% Determine whether the following coordinate and velocity field           %
% structures are correct for MATfluids.                                   %
%                                                                         %
% >> clear coord vel;                                                     %
% >> coord.x = (-2:0.05:2).';                                             %
% >> coord.y = (0:0.1:1).';                                               %
% >> coord.t = (0:0.01:1).';                                              %
% >> vel.w   = rand(length(coord.x), 5, length(coord.t));                 %
% >> examineFlowData(coord, vel);                                         %
% Error using examineFlowData (line 224)                                  %
% The fields of the velocity field structure have an incompatible size    %
% with the provided coordinate structure.                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function examineFlowData(coord, vel)


%% PARSE INPUTS

% Input defaults.
% N/A

% Input checks.
% N/A

% Parse the inputs.
% N/A

% Additional verifications.
nargoutchk(0,0);
examineCoord(coord);
examineVel(vel);


%% EXAMINE COORDINATE AND VELOCITY FIELD COMBINATION

% Get the fields in the coordinate and velocity field structures.
fieldsCoord    = fieldnames(coord);
numFieldsCoord = length(fieldsCoord);
fieldsVel      = fieldnames(vel);
numFieldsVel   = length(fieldsVel);

% Determine the sizes of the coordinate fields.
sizeCoord = zeros(1, numFieldsCoord);
for k = 1:numFieldsCoord
    sizeCoord(k) = length(coord.(fieldsCoord{k}));
end
sizeCoord = sizeCoord(sizeCoord > 0);
sizeCoordNoUnit = sizeCoord(sizeCoord > 1);
clear fieldsCoord k numFieldsCoord;

% Return an error if the size of the velocity field is incompatible with
% the provided coordinates.
C = zeros(3,numFieldsVel);
for k = 1:numFieldsVel
    C(1,k) = isequal(sizeCoord, size(vel.(fieldsVel{k})))               ...
          || isequal(sizeCoord, [size(vel.(fieldsVel{k})) 1]);
    C(2,k) = isequal(sizeCoordNoUnit, size(vel.(fieldsVel{k})));
    C(3,k) = isempty(vel.(fieldsVel{k}))                                ...
          || isscalar(vel.(fieldsVel{k}));
end
if ~min(C(1,:) | C(2,:) | C(3,:))
    error('flowdata:incompatibleSize',                                  ...
         ['The fields of the velocity field structure have an incompat' ...
          'ible size with the provided coordinate structure.']);
end
clear C fieldsVel numFieldsVel sizeCoord sizeCoordNoUnit;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  NOTES                                  %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Line(s) 209                                                             %
% * Empty coordinate fields are permitted by MATfluids.                   %
%                                                                         %
% Line(s) 213-228                                                         %
% * The condition C(1) takes care of the case where one of the fields in  %
%   the coordinate structure has unit size and the unit dimension is also %
%   contained in the velocity field. If the unit size corresponds to the  %
%   last field (typically t or z), then a one should be added to the size %
%   of the velocity field for proper comparison. Alternatively, the final %
%   element in sizeCoord could be ignored.                                %
% * The condition C(2) takes care of the case where one of the fields in  %
%   the coordinate structure has unit size and the unit dimension is not  %
%   contained in the velocity field. In this case, all unit sized fields  %
%   in sizeCoord should be ignored for the comparison.                    %
% * The condition C(3) allows a velocity component to be specified as an  %
%   empty array or a scalar.                                              %
% * Note that only one field in the velocity field structure need be      %
%   evaluated since the examineVel function ensures that all the velocity %
%   field components have the same array size.                            %
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
% 2022/02/23 -- (GDL) Removed message suppression in file, prefer line.   %
% 2022/02/22 -- (GDL) Moved change log and future updates to bottom,      %
%                     reformatted notes.                                  %
% 2021/06/05 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
