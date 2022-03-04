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
% val = examineCoord(coord);                                              %
% [val, msg] = examineCoord(coord);                                       %
% [___] = examineCoord(coord, fields);                                    %
% [___] = examineCoord(coord, fields, 'all');                             %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Examine whether a coordinate structure is properly formatted for        %
% MATfluids. The function will output true if the coordinate structure is %
% properly formatted or false otherwise. By default, the permitted fields %
% in the coordinate structure can be any subset of ["x" "y" "z" "t"] as   %
% long as the order is respected. The user can instead choose to restrict %
% the permitted field names to a smaller subset (ex., only ["x" "t"]). An %
% optional output message can be used to describe what is wrong with the  %
% coordinate structure.                                                   %
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
% 'fields'       STRING ARRAY (1 X N)                                     %
%                Default: ["x" "y" "z" "t"]                               %
%              ~ Permitted field names. The fields in the coordinate      %
%                structure are compared against the permitted fields. The %
%                string array may contain any subset of the fields ["x"   %
%                "y" "z" "t"]. Field names cannot be repeated. The        %
%                ordering of the listed fields is corrected for by the    %
%                function.                                                %
% ----------------------------------------------------------------------- %
% 'all'          CASE-INSENSITIVE CHARACTER ARRAY                         %
%              ~ Specify 'all' to determine whether all the specified     %
%                field names are present in the coordinate structure.     %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'val'          LOGICAL SCALAR                                           %
%              ~ Truth value. The value will return true (1) if the       %
%                coordinate structure is identified as valid and false    %
%                (0) otherwise.                                           %
% ----------------------------------------------------------------------- %
% 'msg'          STRING SCALAR                                            %
%              ~ Coordinate structure message. The string will return     %
%                "valid" if the coordinate structure is identified as     %
%                valid. Otherwise, the string will display a message      %
%                describing an error of the coordinate structure.         %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Determine whether the elements of the following coordinate structure is %
% correct for MATfluids.                                                  %
%                                                                         %
% >> coord.x = (-10:0.1:10).';                                            %
% >> coord.y = (0:0.1:5).';                                               %
% >> coord.t = (0:0.1:2).';                                               %
% >> val     = examineCoord(coord);                                       %
% >> disp(val);                                                           %
%    1                                                                    %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Determine whether the elements of the following coordinate structure is %
% correct for MATfluids.                                                  %
%                                                                         %
% >> coord.t    = (0:0.1:2).';                                            %
% >> coord.y    = (0:0.1:5).';                                            %
% >> coord.x    = (-10:0.1:10).';                                         %
% >> [val, msg] = examineCoord(coord);                                    %
% >> disp(val);                                                           %
%    0                                                                    %
% >> disp(msg);                                                           %
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
% >> coord.x    = [1 2 4 6 5].';                                          %
% >> coord.y    = (0:0.1:5).';                                            %
% >> coord.t    = (0:0.1:2).';                                            %
% >> [val, msg] = examineCoord(coord);                                    %
% >> disp(val);                                                           %
%    0                                                                    %
% >> disp(msg);                                                           %
% At least one field in the coordinate structure does not contain a       %
% strictly monotonically increasing/decreasing real-valued column vector. %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 4                                                               %
%                                                                         %
% Determine whether the elements of the following coordinate structure is %
% correct for MATfluids. Restrict the permitted field names to ["x" "y"   %
% "t"].                                                                   %
%                                                                         %
% >> coord.x    = (-10:0.1:10).';                                         %
% >> coord.y    = (0:0.1:5).';                                            %
% >> coord.z    = (0:0.1:20).';                                           %
% >> coord.t    = (0:0.1:2).';                                            %
% >> [val, msg] = examineCoord(coord, ["x" "y" "t"]);                     %
% >> disp(val);                                                           %
%    0                                                                    %
% >> disp(msg);                                                           %
% At least one field name is not valid in the coordinate structure.       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [val, varargout] = examineCoord(coord, varargin)


%% PARSE INPUTS

% Input defaults.
default.fields = ["x" "y" "z" "t"];
default.all    = 0;

% Input checks.
check.coord  = @(x) validateattributes(x,                               ...
                    {'struct'},                                         ...
                    {'size', [1 1]});
check.fields = @(x) validateattributes(x,                               ...
                    {'string'},                                         ...
                    {'nonempty', 'row'});
check.all    = @(x) any(validatestring(x, {'all'}));

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser , 'coord'  ,                  check.coord  );
addOptional ( hParser , 'fields' , default.fields , check.fields );
addOptional ( hParser , 'all'    , default.all    , check.all    );
parse(hParser, coord, varargin{:});
clear check;

% Additional verifications.
nargoutchk(0,2);

% Verify whether the search fields are valid.
validFields = hParser.Results.fields;
[Lia, Locb] = ismember(validFields, default.fields);
if ~min(Lia)
    error('fields:invalidFields',                                       ...
         ['At least one of the specified search fields is invalid. Onl' ...
          'y a subset of (x,y,z,t) may be used.']);
end
if length(Locb) ~= length(unique(Locb))
    error('fields:repeatedFields',                                      ...
         ['At least one of the specified search fields appears more th' ...
          'an once.']);
end

% Verify, or ensure, that the search fields appear in the proper order.
validFields = default.fields(sort(Locb));
clear default;


%% EXAMINE COORDINATE STRUCTURE

% Get the fields in the coordinate structure.
coordFields = fieldnames(coord);
numFields   = length(coordFields);

% Examine the validity of the fields in the coordinate structure.
[Lia, Locb] = ismember(coordFields, validFields);

% Return false if one of the fields is invalid.
if ~min(Lia)
    val = false;
    msg = ['At least one field name is not valid in the coordinate str' ...
           'ucture.'];
    varargout{1} = sprintf('%s', msg);
    clear msg;
    return;
end
clear Lia;

% Return false if one of the fields is not in the proper order or of one of
% the fields is repeated.
if min(diff(Locb)) <= 0
    val = false;
    msg = ['At least one field in the coordinate structure is not in t' ...
           'he proper order.'];
    varargout{1} = sprintf('%s', msg);
    clear msg;
    return;
end
clear Locb;

% When the 'all' option is checked, return false if all the fields are not
% present in the coordinate structure.
if hParser.Results.all
    if length(coordFields) ~= length(validFields)
        val = false;
        msg = ['At least one field is missing from the coordinate stru' ...
               'cture.'];
        varargout{1} = sprintf('%s', msg);
        clear msg;
        return;
    end
end
clear validFields;

% Return false if one of the coordinate fields is not a column vector.
C = zeros(3,numFields);
for k = 1:numFields
    C(1,k) = isrealnum(coord.(coordFields{k}), 'all');
    C(2,k) = iscolumn(coord.(coordFields{k}));
    C(3,k) = isempty(coord.(coordFields{k}))                            ...
          || isscalar(coord.(coordFields{k}));
end
if ~min((C(1,:) & C(2,:)) | C(3,:))
    val = false;
    msg = ['At least one field in the coordinate structure does not co' ...
           'ntain a strictly monotonically increasing/decreasing real-' ...
           'valued column vector.'];
    varargout{1} = sprintf('%s', msg);
    clear msg;
    return;
end

% Return false if one of the coordinate fields is not strictly
% monotonically increasing/decreasing.
for k = 1:numFields
    tmp  = diff(coord.(coordFields{k}));
    C(1,k) = isposreal(tmp, 'all');
    C(2,k) = isnegreal(tmp, 'all');
end
if ~min(C(1,:) | C(2,:) | C(3,:))
    val = false;
    msg = ['At least one field in the coordinate structure does not co' ...
           'ntain a strictly monotonically increasing/decreasing real-' ...
           'valued column vector.'];
    varargout{1} = sprintf('%s', msg);
    return;
end
clear C coordFields k numFields tmp;

% If the code reached this point, all tests have passed.
val          = true;
varargout{1} = "valid";


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
% 2022/03/04 -- (GDL) Added 'all' option.                                 %
% 2022/03/04 -- (GDL) Changed function behaviour.                         %
% 2022/02/23 -- (GDL) Removed message suppression in file, prefer line.   %
% 2022/02/22 -- (GDL) Moved change log and future updates to bottom,      %
%                     reformatted notes.                                  %
% 2021/06/06 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
