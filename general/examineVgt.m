%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% examineVgt                                                              %
% MATfluids Structure Interrogation                                       %
% Examine contents of velocity gradient structure                         %
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
% val = examineVgt(vgt);                                                  %
% [val, msg] = examineVgt(vgt);                                           %
% [___] = examineVgt(vgt, fields);                                        %
% [___] = examineVgt(vgt, fields, 'exact');                               %
% [___] = examineVgt(___, Name, Value);                                   %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Examine whether a velocity gradient structure is properly formatted for %
% MATfluids. The function will output true if the structure is properly   %
% formatted or false otherwise. By default, the permitted fields in the   %
% velocity gradient structure can be any subset of ["ux" "uy" "uz" "vx"   %
% "vy" "vz" "wx" "wy" "wz"] as long as the order is respected. The user   %
% can instead choose to restrict the permitted field names to a smaller   %
% subset (ex., only ["ux" "uy" "vx" "vy"]). An optional output message    %
% can be used to describe what is wrong with the velocity gradient        %
% structure.                                                              %
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
% 'vgt'          STRUCT ARRAY (1 X 1)                                     %
%              ~ Velocity gradient tensor. The structure array may        %
%                contain the following fields:                            %
%                1) ux : Gradient of u component in the x direction       %
%                2) uy : Gradient of u component in the y direction       %
%                3) uz : Gradient of u component in the z direction       %
%                4) vx : Gradient of v component in the x direction       %
%                5) vy : Gradient of v component in the y direction       %
%                6) vz : Gradient of v component in the z direction       %
%                7) wx : Gradient of w component in the x direction       %
%                8) wy : Gradient of w component in the y direction       %
%                9) wz : Gradient of w component in the z direction       %
%                Each field is a one to four dimensional array, in ndgrid %
%                format (i.e., dimensions represent [x y z t]). All       %
%                constituent arrays must have the same size. Not all      %
%                fields need be present in the velocity gradient          %
%                structure, however the order must follow [ux uy uz vx vy %
%                vz wx wy wz]. For example, the velocity gradient         %
%                structure may contain fields ordered as [ux uy vx vy]    %
%                but not as [vx vy ux uy]. Empty, scalar and complex      %
%                arrays are permitted.                                    %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'fields'       STRING ARRAY (1 X N)                                     %
%                Default: ["ux" "uy" "uz" "vx" "vy" "vz" "wx" "wy" "wz"]  %
%              ~ Permitted field names. The fields in the velocity        %
%                gradient structure are compared against the permitted    %
%                fields. The string array may contain any subset of the   %
%                fields ["ux" "uy" "uz" "vx" "vy" "vz" "wx" "wy" "wz"].   %
%                Field names cannot be repeated. The ordering of the      %
%                listed fields is corrected for by the function.          %
% ----------------------------------------------------------------------- %
% 'exact'        CASE-INSENSITIVE CHARACTER ARRAY                         %
%              ~ Specify 'exact' to determine whether all and only the    %
%                specified field names are present in the velocity        %
%                gradient structure.                                      %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% 'scalar'       CASE-INSENSITIVE CHARACTER ARRAY                         %
%                Default: 'on'                                            %
%              ~ Specify whether the velocity gradient components can     %
%                have scalar or empty values. Specifying 'on' will allow  %
%                it (default) and specifying 'off' will disallow it.      %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'val'          LOGICAL SCALAR                                           %
%              ~ Truth value. The value will return true (1) if the       %
%                velocity gradient structure is identified as valid and   %
%                false (0) otherwise.                                     %
% ----------------------------------------------------------------------- %
% 'msg'          STRING SCALAR                                            %
%              ~ Velocity gradient structure message. The string will     %
%                return "valid" if the velocity gradient structure is     %
%                identified as valid. Otherwise, the string will display  %
%                a message describing an error of the velocity gradient   %
%                structure.                                               %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Determine whether the elements of the following velocity gradient       %
% structure is correct for MATfluids.                                     %
%                                                                         %
% >> vgt.ux = 1;                                                          %
% >> vgt.uy = [1 2 3; 4 5 6; 7 8 9];                                      %
% >> vgt.vx = [1 2 1; 2 1 2; 1 2 1];                                      %
% >> vgt.vy = [];                                                         %
% >> val = examineVgt(vgt);                                               %
% >> disp(val);                                                           %
% >> 1                                                                    %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Determine whether the elements of the following velocity gradient       %
% structure is correct for MATfluids.                                     %
%                                                                         %
% >> vgt.ux = [-1 1; -2 2];                                               %
% >> vgt.uy = [1 2 3; 4 5 6; 7 8 9];                                      %
% >> vgt.vx = [1 2 1; 2 1 2; 1 2 1];                                      %
% >> vgt.vy = [];                                                         %
% >> [val, msg] = examineVgt(vgt);                                        %
% >> disp(val);                                                           %
%    0                                                                    %
% >> disp(msg);                                                           %
% At least one field in the velocity gradient structure has a different   %
% array size than the others.                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [val, varargout] = examineVgt(vgt, varargin)


%% PARSE INPUTS

% Input defaults.
default.fields = ["ux" "uy" "uz" "vx" "vy" "vz" "wx" "wy" "wz"];
default.exact  = 0;
default.scalar = 'on';

% Input checks.
check.vgt    = @(x) validateattributes(x,                               ...
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
addRequired ( hParser , 'vgt'    ,                  check.vgt    );
addOptional ( hParser , 'fields' , default.fields , check.fields );
addOptional ( hParser , 'exact'  , default.exact  , check.exact  );
addParameter( hParser , 'scalar' , default.scalar , check.scalar );
parse(hParser, vgt, varargin{:});
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
          'y a subset of (ux,uy,uz,vx,vy,vz,wx,wy,wz) may be used.']);
end
if length(Locb) ~= length(unique(Locb))
    error('fields:repeatedFields',                                      ...
         ['At least one of the specified search fields appears more th' ...
          'an once.']);
end

% Verify, or ensure, that the search fields appear in the proper order.
validFields = default.fields(sort(Locb));
clear default;

%% EXAMINE VELOCITY GRADIENT STRUCTURE

% Get the fields in the velocity gradient structure.
vgtFields = fieldnames(vgt);
numFields = length(vgtFields);

% Examine the validity of the fields in the velocity gradient structure.
[Lia, Locb] = ismember(vgtFields, validFields);

% Return false if one of the fields is invalid.
if ~min(Lia)
    val = false;
    msg = ['At least one field name is not valid in the velocity gradi' ...
           'ent structure.'];
    varargout{1} = sprintf('%s', msg);
    clear msg;
    return;
end
clear Lia;

% Return false if one of the fields is not in the proper order or if one of
% the fields is repeated.
if (min(diff(Locb)) <= 0) || (length(Locb) ~= length(unique(Locb)))
    val = false;
    msg = ['At least one field in the velocity gradient structure is n' ...
           'ot in the proper order or appears more than once.'];
    varargout{1} = sprintf('%s', msg);
    clear msg;
    return;
end
clear Locb;

% When the 'exact' option is checked, return false if the specified fields
% do not exactly match the velocity gradient structure.
if hParser.Results.exact
    if length(vgtFields) ~= length(validFields)
        val = false;
        msg = ['At least one field is missing from the velocity gradie' ...
               'nt structure.'];
        varargout{1} = sprintf('%s', msg);
        clear msg;
        return;
    end
end
clear validFields;

% Return false if one of the fields is not a numeric array.
C = zeros(2,numFields);
for k = 1:numFields
    C(1,k) = islogiconum(vgt.(vgtFields{k}));
    C(2,k) = isempty(vgt.(vgtFields{k})) || isscalar(vgt.(vgtFields{k}));
end
if ~min(C(1,:) | C(2,:))
    val = false;
    msg = ['At least one field in the velocity gradient structure does' ...
          ' not contain a numeric array.'];
    varargout{1} = sprintf('%s', msg);
    clear msg;
    return;
end

% Return an error if one of the fields has a different size.
flag = 1;
allowScalar = strcmpi(hParser.Results.scalar, 'on');
for k = 1:numFields
    % Continue if the velocity gradient component is empty or scalar.
    if C(2,k) && allowScalar, continue; end
    
    % Obtain the array size of the velocity gradient component for the
    % first time.
    if flag
        sz = size(vgt.(vgtFields{k}));
        flag = 0;
        continue;
    end
    
    % Compare array sizes for velocity gradient components.
    if ~isequal(sz, size(vgt.(vgtFields{k})))
        val = false;
        msg = ['At least one field in the velocity gradient structure ' ...
               'has a different array size than the others.'];
        varargout{1} = sprintf('%s', msg);
        clear msg;
        return;
    end
end
clear allowScalar C flag k numFields sz vgtFields;

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
% 2023/09/22 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
