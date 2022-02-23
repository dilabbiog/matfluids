%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% isthru                                                                  %
% Data Interrogation                                                      %
% Examine data for a type through cells and structs.                      %
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
% val = isthru(S, @fun);                                                  %
% val = isthru(___, Option, Name, Value);                                 %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
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
% N/A                                                                     %
%                                                                         %
% ======================================================================= %
% Input Arguments (Required):                                             %
% ----------------------------------------------------------------------- %
% 'fun'          FUNCTION HANDLE                                          %
%              ~ Function handle of data type interrogator. For example,  %
%                this could be @isint, @isnumeric, @isposreal, @isstring, %
%                @isyes, etc.                                             %
% ----------------------------------------------------------------------- %
% 'S'            N-DIMENSIONAL CELL/STRUCT/NUMERIC ARRAY                  %
%              ~ Input array. The elements of this array will be examined %
%                for genuine integers. Data types that are not logical or %
%                numeric are accepted, but will return zero.              %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% N/A            N/A                                                      %
%              ~ The valid options depend on the data type interrogator   %
%                specified by 'fun'.                                      %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% N/A            N/A                                                      %
%              ~ The valid name-value pairs depend on the data type       %
%                interrogator specified by 'fun'.                         %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'val'          LOGICAL N-DIMENSIONAL CELL/STRUCT/NUMERIC ARRAY/SCALAR   %
%              ~ Output array. The elements of this array address whether %
%                the elements of the input array are of the specified     %
%                data type (1) or not (0).                                %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Create a variable x that consists of nested cell arrays and structs     %
% whose elements are ultimately arrays of some sort (character, numeric,  %
% etc.). Determine whether the elements are positive real.                %
%                                                                         %
% >> x      = cell(3,2);                                                  %
% >> x{1,1} = struct('a', [1 2; -1 -2], 'b', ['a' 'b' 'c']);              %
% >> x{2,1} = [-2.22 -Inf NaN; Inf 10^-7 2];                              %
% >> x{3,1} = ["Hello" "Goodbye"]';                                       %
% >> x{1,2} = struct('d', {[5 4 3; -2 1 0]; pi/2}, 'e', {-2; 20});        %
% >> x{2,2} = "Another cell!";                                            %
% >> x{3,2} = -0.01;                                                      %
% >> val    = isthru(x, @isposreal);                                      %
% >> disp(val{1,1}.a);                                                    %
%    1   1                                                                %
%    0   0                                                                %
% >> disp(val{1,1}.b);                                                    %
%      0                                                                  %
% >> disp(val{1,2}(1,1).d);                                               %
%    1   1   1                                                            %
%    0   1   0                                                            %
% >> disp(val{1,2}(2,1).d);                                               %
%    1                                                                    %
% >> disp(val{1,2}(1,1).e);                                               %
%    0                                                                    %
% >> disp(val{1,2}(2,1).e);                                               %
%    1                                                                    %
% >> disp(val{2,1});                                                      %
%    0   0   1                                                            %
%    1   1   1                                                            %
% >> disp(val{2,2});                                                      %
%      0                                                                  %
% >> disp(val{3,1});                                                      %
%      0                                                                  %
% >> disp(val{3,2});                                                      %
%      0                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [val] = isthru(S, fun, varargin)


%% PARSE INPUTS

% Input defaults.
% N/A

% Input checks.
check.fun = @(x) validateattributes(x,                                  ...
                 {'function_handle'},                                   ...
                 {});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'S'               );
addRequired ( hParser, 'fun' , check.fun );
parse(hParser, S, fun, varargin{:});
clear check default;

% Additional verifications.
nargoutchk(0,1);


%% DATA TYPE INTERROGATION

% Determine the variable type of 'S'.
if ischar(S) || islogical(S) || isnumeric(S) || isstring(S)
    type = 'array';
elseif isstruct(S)
    type = 'struct';
elseif iscell(S)
    type = 'cell';
end

% Apply the interrogation function.
switch type
    case 'array'
        val = fun(S, varargin{:});
    case 'cell'
        val = cell(size(S));
        for k = 1:numel(S)
            val{ind2sub(size(S),k)}                                     ...
                = isthru(S{ind2sub(size(S),k)}, fun, varargin{:});
        end
    case 'struct'
        fields = fieldnames(S);
        val = cell2struct(cell(length(fields),1),fields);
        val = repmat(val, size(S));
        for k1 = 1:numel(S)
            for k2 = 1:length(fields)
                val(ind2sub(size(S),k1)).(fields{k2})                   ...
                    = isthru(S(ind2sub(size(S),k1)).(fields{k2}), fun,  ...
                             varargin{:});
            end
        end
end


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
%#ok<*N/A>                                                                %
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
% 2022/02/22 -- (GDL) Moved change log and future updates to bottom,      %
%                     reformatted notes.                                  %
% 2021/06/04 -- (GDL) Added nargoutchk and future updates comments.       %
% 2021/05/27 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
