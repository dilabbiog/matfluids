%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% examineVel                                                              %
% MATfluids Structure Interrogation                                       %
% Examine contents of velocity field structure                            %
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
% 2021/06/05 -- (GDL) Beta version of the code finalized.                 %
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
% examineVel(vel);                                                        %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Examine whether a velocity field structure is properly formatted for    %
% MATfluids. The function has no specific output, it will simply throw an %
% error if the velocity field structure is not correctly formatted.       %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2019b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
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
% Determine whether the elements of the following velocity field          %
% structure is correct for MATfluids.                                     %
%                                                                         %
% >> clear vel;                                                           %
% >> vel.u = 1;                                                           %
% >> vel.v = [1 2 3; 4 5 6; 7 8 9];                                       %
% >> vel.w = [1 2 1; 2 1 2; 1 2 1];                                       %
% >> examineVel(vel);                                                     %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% >> clear vel;                                                           %
% >> vel.u = [-1 1; -2 2];                                                %
% >> vel.v = [1 2 3; 4 5 6; 7 8 9];                                       %
% >> vel.w = [];                                                          %
% >> examineVel(vel);                                                     %
% Error using examineVel (line 242)                                       %
% At least one field in the velocity field structure has a different      %
% array size than the others.                                             %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 3                                                               %
%                                                                         %
% >> clear vel;                                                           %
% >> vel.u = [-1 1; -2 2];                                                %
% >> vel.v = [1 2; 3 4];                                                  %
% >> vel.w = [];                                                          %
% >> examineVel(vel);                                                     %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 4                                                               %
%                                                                         %
% >> clear vel;                                                           %
% >> vel.v = [-1 1; -2 2];                                                %
% >> vel.u = [1 2; 3 4];                                                  %
% >> examineVel(vel);                                                     %
% Error using examineVel (line 206)                                       %
% At least one field in the velocity field structure is not in the proper %
% order.                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function examineVel(vel)


%% PARSE INPUTS

% Input defaults.
% N/A

% Input checks.
check.vel = @(x) validateattributes(x,                                  ...
                 {'struct'},                                            ...
                 {'size', [1 1]});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'vel' , check.vel );
parse(hParser, vel);
clear check;

% Additional verifications.
nargoutchk(0,0);


%% EXAMINE VELOCITY FIELD STRUCTURE

% List of valid fields.
validFields = ["u"; "v"; "w"];

% Get the fields in the velocity field structure.
fields    = fieldnames(vel);
numFields = length(fields);

% Examine the validity of the fields in the velocity field structure.
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
   error('vel:invalidFields',                                           ...
        ['At least one field name is not valid in the velocity field s' ...
         'tructure.']);
end
clear val;

% Return an error if one of the fields is not in the proper order.
if min(diff(idx)) <= 0
   error('vel:unorderedFields',                                         ...
        ['At least one field in the velocity field structure is not in' ...
         ' the proper order.']);
end
clear idx;

% Return an error if one of the fields is not a real-valued array.
C = zeros(2,numFields);
for k = 1:numFields
    C(1,k) = isrealnum(vel.(fields{k}), 'all');
    C(2,k) = isempty(vel.(fields{k}));
end
if ~min(C(1,:) | C(2,:))
    error('vel:notReal',                                                ...
         ['At least one field in the velocity field structure does not' ...
          'contain a real-valued array.']);
end

% Return an error if one of the fields has a different size.
flag = 1;
for k = 1:numFields
    % Continue if the velocity component is empty or scalar.
    C(1,k) = isscalar(vel.(fields{k}));
    if C(1,k) | C(2,k), continue; end
    
    % Obtain the array size of the velocity field component for the first
    % time.
    if flag
        tmp  = size(vel.(fields{k}));
        flag = 0;
        continue;
    end
    
    % Compare array sizes for nonempty and nonscalar velocity field
    % components.
    if ~isequal(tmp, size(vel.(fields{k})))
        error('vel:invalidSize',                                        ...
             ['At least one field in the velocity field structure has ' ...
              'a different array size than the others.']);
    end
end
clear C fields flag k numFields tmp;


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
