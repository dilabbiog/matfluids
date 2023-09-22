%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% examineVel                                                              %
% MATfluids Structure Interrogation                                       %
% Examine contents of velocity structure                                  %
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
% val = examineVel(vel);                                                  %
% [val, msg] = examineVel(vel);                                           %
% [___] = examineVel(vel, fields);                                        %
% [___] = examineVel(vel, fields, 'exact');                               %
% [___] = examineVel(___, Name, Value);                                   %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Examine whether a velocity structure is properly formatted for          %
% MATfluids. The function will output true if the structure is properly   %
% formatted or false otherwise. By default, the permitted fields in the   %
% velocity structure can be any subset of ["u" "v" "w"] as long as the    %
% order is respected. The user can instead choose to restrict the         %
% permitted field names to a smaller subset (ex., only ["u" "v"]). An     %
% optional output message can be used to describe what is wrong with the  %
% velocity structure.                                                     %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2019b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
% islogiconum                                                             %
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
% 'vel'          STRUCT ARRAY (1 X 1)                                     %
%              ~ Velocity field. The structure array may contain the      %
%                following fields:                                        %
%                1) u : Velocity component in the x direction             %
%                2) v : Velocity component in the y direction             %
%                3) w : Velocity component in the z direction             %
%                Each field is a one to four dimensional array, in ndgrid %
%                format (i.e., dimensions represent [x y z t]). All       %
%                constituent arrays must have the same size. Not all      %
%                fields need be present in the velocity structure,        %
%                however the order must follow [u v w]. For example, the  %
%                velocity structure may contain fields ordered as [u w]   %
%                but not as [w u]. Empty, scalar and complex arrays are   %
%                permitted.                                               %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'fields'       STRING ARRAY (1 X N)                                     %
%                Default: ["u" "v" "w"]                                   %
%              ~ Permitted field names. The fields in the velocity        %
%                structure are compared against the permitted fields. The %
%                string array may contain any subset of the fields ["u"   %
%                "v" "w"]. Field names cannot be repeated. The ordering   %
%                of the listed fields is corrected for by the function.   %
% ----------------------------------------------------------------------- %
% 'exact'        CASE-INSENSITIVE CHARACTER ARRAY                         %
%              ~ Specify 'exact' to determine whether all and only the    %
%                specified field names are present in the velocity        %
%                structure.                                               %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% 'scalar'       CASE-INSENSITIVE CHARACTER ARRAY                         %
%                Default: 'on'                                            %
%              ~ Specify whether the velocity components can have scalar  %
%                or empty values. Specifying 'on' will allow it (default) %
%                and specifying 'off' will disallow it.                   %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'val'          LOGICAL SCALAR                                           %
%              ~ Truth value. The value will return true (1) if the       %
%                velocity structure is identified as valid and false (0)  %
%                otherwise.                                               %
% ----------------------------------------------------------------------- %
% 'msg'          STRING SCALAR                                            %
%              ~ Velocity structure message. The string will return       %
%                "valid" if the velocity structure is identified as       %
%                valid. Otherwise, the string will display a message      %
%                describing an error of the velocity structure.           %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Determine whether the elements of the following velocity structure is   %
% correct for MATfluids.                                                  %
%                                                                         %
% >> vel.u = 1;                                                           %
% >> vel.v = [1 2 3; 4 5 6; 7 8 9];                                       %
% >> vel.w = [1 2 1; 2 1 2; 1 2 1];                                       %
% >> val = examineVel(vel);                                               %
% >> disp(val);                                                           %
% >> 1                                                                    %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Determine whether the elements of the following velocity structure is   %
% correct for MATfluids.                                                  %
%                                                                         %
% >> vel.u = [-1 1; -2 2];                                                %
% >> vel.v = [1 2 3; 4 5 6; 7 8 9];                                       %
% >> vel.w = [];                                                          %
% >> [val, msg] = examineVel(vel);                                        %
% >> disp(val);                                                           %
%    0                                                                    %
% >> disp(msg);                                                           %
% At least one field in the velocity structure has a different array size %
% than the others.                                                        %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 3                                                               %
%                                                                         %
% Determine whether the elements of the following velocity structure is   %
% correct for MATfluids.                                                  %
%                                                                         %
% >> vel.u = [-1 1; -2 2];                                                %
% >> vel.v = [1 2; 3 4];                                                  %
% >> vel.w = [];                                                          %
% >> [val, msg] = examineVel(vel);                                        %
% >> disp(val);                                                           %
%    1                                                                    %
% >> disp(msg);                                                           %
% valid                                                                   %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 4                                                               %
%                                                                         %
% Determine whether the elements of the following velocity structure is   %
% correct for MATfluids.                                                  %
%                                                                         %
% >> vel.v = [-1 1; -2 2];                                                %
% >> vel.u = [1 2; 3 4];                                                  %
% >> [val, msg] = examineVel(vel);                                        %
% >> disp(val);                                                           %
%    0                                                                    %
% >> disp(msg);                                                           %
% At least one field in the velocity structure is not in the proper order %
% or appears more than once.                                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [val, varargout] = examineVel(vel, varargin)


%% PARSE INPUTS

% Input defaults.
default.fields = ["u" "v" "w"];
default.exact  = 0;
default.scalar = 'on';

% Input checks.
check.vel    = @(x) validateattributes(x,                               ...
                    {'struct'},                                         ...
                    {'size', [1 1]});
check.fields = @(x) validateattributes(x,                               ...
                    {'string'},                                         ...
                    {'nonempty', 'row'});
check.exact  = @(x) any(validatestring(x, {'exact'}));
check.scalar = @(x) any(validatestring(x,                               ...
                       {'on' 'off'}));

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser , 'vel'    ,                  check.vel    );
addOptional ( hParser , 'fields' , default.fields , check.fields );
addOptional ( hParser , 'exact'  , default.exact  , check.exact  );
addParameter( hParser , 'scalar' , default.scalar , check.scalar );
parse(hParser, vel, varargin{:});
clear check;

% Additional verifications.
narginchk(1,5);
nargoutchk(0,2);

% Verify whether the search fields are valid.
validFields = hParser.Results.fields;
[Lia, Locb] = ismember(validFields, default.fields);
if ~min(Lia)
    error('fields:invalidFields',                                       ...
         ['At least one of the specified search fields is invalid. Onl' ...
          'y a subset of (u,v,w) may be used.']);
end
if length(Locb) ~= length(unique(Locb))
    error('fields:repeatedFields',                                      ...
         ['At least one of the specified search fields appears more th' ...
          'an once.']);
end

% Verify, or ensure, that the search fields appear in the proper order.
validFields = default.fields(sort(Locb));
clear default;

%% EXAMINE VELOCITY STRUCTURE

% Get the fields in the velocity structure.
velFields = fieldnames(vel);
numFields = length(velFields);

% Examine the validity of the fields in the velocity structure.
[Lia, Locb] = ismember(velFields, validFields);

% Return false if one of the fields is invalid.
if ~min(Lia)
    val = false;
    msg = ['At least one field name is not valid in the velocity struc' ...
           'ture.'];
    varargout{1} = sprintf('%s', msg);
    clear msg;
    return;
end
clear Lia;

% Return false if one of the fields is not in the proper order or if one of
% the fields is repeated.
if (min(diff(Locb)) <= 0) || (length(Locb) ~= length(unique(Locb)))
    val = false;
    msg = ['At least one field in the velocity structure is not in the' ...
           ' proper order or appears more than once.'];
    varargout{1} = sprintf('%s', msg);
    clear msg;
    return;
end
clear Locb;

% When the 'exact' option is checked, return false if the specified fields
% do not exactly match the velocity structure.
if hParser.Results.exact
    if length(velFields) ~= length(validFields)
        val = false;
        msg = ['At least one field is missing from the velocity struct' ...
               'ure.'];
        varargout{1} = sprintf('%s', msg);
        clear msg;
        return;
    end
end
clear validFields;

% Return false if one of the fields is not a numeric array.
C = zeros(2,numFields);
for k = 1:numFields
    C(1,k) = islogiconum(vel.(velFields{k}));
    C(2,k) = isempty(vel.(velFields{k})) || isscalar(vel.(velFields{k}));
end
if ~min(C(1,:) | C(2,:))
    val = false;
    msg = ['At least one field in the velocity structure does not cont' ...
          'ain a numeric array.'];
    varargout{1} = sprintf('%s', msg);
    clear msg;
    return;
end

% Return an error if one of the fields has a different size.
flag = 1;
allowScalar = strcmpi(hParser.Results.scalar, 'on');
for k = 1:numFields
    % Continue if the velocity component is empty or scalar.
    if C(2,k) && allowScalar, continue; end
    
    % Obtain the array size of the velocity component for the first time.
    if flag
        sz = size(vel.(velFields{k}));
        flag = 0;
        continue;
    end
    
    % Compare array sizes for velocity components.
    if ~isequal(sz, size(vel.(velFields{k})))
        val = false;
        msg = ['At least one field in the velocity structure has a dif' ...
               'ferent array size than the others.'];
        varargout{1} = sprintf('%s', msg);
        clear msg;
        return;
    end
end
clear allowScalar C flag k numFields sz velFields;

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
% 2023/09/22 -- (GDL) Updated to be more like examineCoord.               %
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
