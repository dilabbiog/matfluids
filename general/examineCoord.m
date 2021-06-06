%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% examineCoord                                                            %
% MATfluids Structure Interrogation                                       %
% Examine contents of coordinate structure                                %
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
% CHANGE LOG                                                              %
%                                                                         %
% 2021/06/06 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Copyright (C) 2021 Giuseppe Di Labbio                                   %
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
% examineCoord(coord);                                                    %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Examine whether a coordinate structure is properly formatted for        %
% MATfluids. The function has no specific output, it will simply throw an %
% error if the coordinate structure is not correctly formatted.           %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2019b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
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
% Determine whether the elements of the following coordinate structure is %
% correct for MATfluids.                                                  %
%                                                                         %
% >> clear coord;                                                         %
% >> coord.x = (-10:0.1:10).';                                            %
% >> coord.y = (0:0.1:5).';                                               %
% >> coord.t = (0:0.1:2).';                                               %
% >> examineCoord(coord);                                                 %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Determine whether the elements of the following coordinate structure is %
% correct for MATfluids.                                                  %
%                                                                         %
% >> clear coord;                                                         %
% >> coord.t = (0:0.1:2).';                                               %
% >> coord.y = (0:0.1:5).';                                               %
% >> coord.x = (-10:0.1:10).';                                            %
% >> examineCoord(coord);                                                 %
% Error using examineCoord (line 205)                                     %
% At least one field in the coordinate structure is not in the proper     %
% order.                                                                  %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 3                                                               %
%                                                                         %
% Determine whether the elements of the following coordinate structure is %
% correct for MATfluids.                                                  %
%                                                                         %
% >> clear coord;                                                         %
% >> coord.x = [1 2 4 6 5].';                                             %
% >> coord.y = (0:0.1:5).';                                               %
% >> coord.t = (0:0.1:2).';                                               %
% >> examineCoord(coord);                                                 %
% Error using examineCoord (line 233)                                     %
% At least one field in the coordinate structure does not contain a       %
% strictly monotonically increasing/decreasing real-valued column vector. %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function examineCoord(coord)


%% PARSE INPUTS

% Input defaults.
% N/A

% Input checks.
check.coord = @(x) validateattributes(x,                                ...
                   {'struct'},                                          ...
                   {'size', [1 1]});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'coord' , check.coord );
parse(hParser, coord);
clear check;

% Additional verifications.
nargoutchk(0,0);


%% EXAMINE COORDINATE STRUCTURE

% List of valid fields.
validFields = ["x"; "y"; "z"; "t"];

% Get the fields in the coordinate structure.
fields    = fieldnames(coord);
numFields = length(fields);

% Examine the validity of the fields in the coordinate structure.
idx = zeros(size(fields));
val = zeros(size(fields));
for k = 1:numFields
    % Determine whether the field is valid and its index in validFields.
    [val(k), idx(k)] = max(strcmpi(fields{k}, validFields));
    
    % If there is no match, ensure the index is set to zero.
    idx(k) = val(k)*idx(k);
end
val = min(val);
clear validFields;

% Return an error if one of the fields is invalid.
if ~val
   error('coord:invalidFields',                                         ...
        ['At least one field name is not valid in the coordinate struc' ...
         'ture.']);
end
clear val;

% Return an error if one of the fields is not in the proper order.
if min(diff(idx)) <= 0
   error('coord:unorderedFields',                                       ...
        ['At least one field in the coordinate structure is not in the' ...
         ' proper order.']);
end
clear idx;

% Return an error if one of the coordinate fields is not a column vector.
C = zeros(3,numFields);
for k = 1:numFields
    C(1,k) = isrealnum(coord.(fields{k}), 'all');
    C(2,k) = iscolumn(coord.(fields{k}));
    C(3,k) = isempty(coord.(fields{k}))                                 ...
          || isscalar(coord.(fields{k}));
end
if ~min((C(1,:) & C(2,:)) | C(3,:))
    error('coord:invalidCoordinates',                                   ...
         ['At least one field in the coordinate structure does not con' ...
          'tain a strictly monotonically increasing/decreasing real-va' ...
          'lued column vector.']);
end

% Return an error if one of the coordinate fields is not strictly
% monotonically increasing/decreasing.
for k = 1:numFields
    tmp  = diff(coord.(fields{k}));
    C(1,k) = isposreal(tmp, 'all');
    C(2,k) = isnegreal(tmp, 'all');
end
if ~min(C(1,:) | C(2,:) | C(3,:))
    error('coord:invalidCoordinates',                                   ...
         ['At least one field in the coordinate structure does not con' ...
          'tain a strictly monotonically increasing/decreasing real-va' ...
          'lued column vector.']);
end
clear C fields k numFields tmp;


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
% * N/A                                                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
