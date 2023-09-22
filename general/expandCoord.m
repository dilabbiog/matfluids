%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% expandCoord                                                             %
% MATfluids Structure Interrogation                                       %
% Expand contents of coordinate structure to most general form            %
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
% Copyright (C) 2023 Giuseppe Di Labbio                                   %
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
% ecoord = expandCoord(coord);                                            %
% ecoord = expandCoord(coord, val);                                       %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Expand a valid MATfluids coordinate structure to include all possible   %
% fields ["x" "y" "z" "t"]. Missing fields are added as empty arrays by   %
% default, though the user can also specify a scalar value.               %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2019b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
% examineCoord                                                            %
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
%                1) x : Spatial coordinate x                              %
%                2) y : Spatial coordinate y                              %
%                3) z : Spatial coordinate z                              %
%                4) t : Temporal coordinate t                             %
%                Each field is a column vector (one-dimensional array)    %
%                that strictly increases/decreases monotonically from the %
%                first row to the last. Not all fields need be present in %
%                the coordinate structure, however the order must follow  %
%                [x y z t]. For example, the coordinate structure may     %
%                contain fields ordered as [x z t] but not as [t x z].    %
%                Both empty and scalar arrays are permitted.              %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'val'          SCALAR                                                   %
%              ~ Input scalar. Missing fields of the coordinate structure %
%                will be assigned this value.                             %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'ecoord'       STRUCT ARRAY (1 X 1)                                     %
%              ~ Expanded coordinates. The structure array contains all   %
%                the following fields:                                    %
%                1) x : Spatial coordinate x                              %
%                2) y : Spatial coordinate y                              %
%                3) z : Spatial coordinate z                              %
%                4) t : Temporal coordinate t                             %
%                Fields available in the input coordinates are copied     %
%                whereas fields that aren't are added as empty arrays.    %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Expand a MATfluids coordinate structure initially containing only the   %
% ["x" "t"] fields.                                                       %
%                                                                         %
% >> coord.x = (0:0.05:10).';                                             %
% >> coord.t = (0:0.5:5).';                                               %
% >> coord   = expandCoord(coord);                                        %
% >> disp(coord);                                                         %
%    x: [201×1 double]                                                    %
%    y: []                                                                %
%    z: []                                                                %
%    t: [11×1 double]                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ecoord] = expandCoord(coord)


%% PARSE INPUTS

% Input defaults.
default.val = [];

% Input checks.
check.coord = @(x) examineCoord(x);
check.val   = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'scalar'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser , 'coord'  ,              check.coord );
addOptional ( hParser , 'val'    , default.val, check.val   );
parse(hParser, coord, varargin{:});
clear check;

% Additional verifications.
narginchk(1,2);
nargoutchk(0,1);


%% INITIALIZATIONS

ecoord = coord;               % Expanded coordinate structure
val    = hParser.Results.val; % Assigned value to missing fields


%% EXPAND COORDINATE STRUCTURE

% Add missing fields.
if ~isfield(ecoord, 'x'), ecoord.x = val; end
if ~isfield(ecoord, 'y'), ecoord.y = val; end
if ~isfield(ecoord, 'z'), ecoord.z = val; end
if ~isfield(ecoord, 't'), ecoord.t = val; end
clear val;

% Order the fields in MATfluids format.
ecoord = orderfields(ecoord, {'x', 'y', 'z', 't'});


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
% 2023/09/22 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
