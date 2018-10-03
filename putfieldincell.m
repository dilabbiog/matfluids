%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%           PUT A NEW STRUCT FIELD INTO A CELL ARRAY OF STRUCTS           %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Giuseppe Di Labbio                                                      %
% Department of Mechanical, Industrial & Aerospace Engineering            %
% Concordia University Montréal, Canada                                   %
%                                                                         %
% Last Update: October 3rd, 2018 by Giuseppe Di Labbio                    %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Copyright (C) 2018 Giuseppe Di Labbio                                   %
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
% C = putfieldincell(data, field, put);                                   %
% C = putfieldincell(data, field, put, index);                            %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function puts a specified field into the structs of a one-         %
% dimensional cell array of structs. The user has the option of selecting %
% where to place the new field in the struct list. By default, the new    %
% field is placed at the end of the field list unless specified otherwise %
% by the user.                                                            %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'C'          - 1D CELL                                                  %
%              - One-dimensional cell array containing the new field.     %
% ----------------------------------------------------------------------- %
% 'data'       - 1D CELL                                                  %
%              - One-dimensional cell array of structs with a field to be %
%                put in.                                                  %
% ----------------------------------------------------------------------- %
% 'field'      - STRING                                                   %
%              - The name of the field to be put in.                      %
% ----------------------------------------------------------------------- %
% 'index'      - REAL SCALAR                                              %
%              - Index in which to place the new field in the list of     %
%                fields.                                                  %
% ----------------------------------------------------------------------- %
% 'put'        - 1D CELL                                                  %
%              - One-dimensional cell array of the field to put in the    %
%                new cell array.                                          %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Given a cell array A of structs with fields 'X', 'Y', and 'Z', put a    %
% new field 'Y2' into the cell array. Place the new field in the third    %
% position in the field lists.                                            %
%                                                                         %
% >> A = cell(5, 1);                                                      %
% >> for k = 1:5                                                          %
% A{k} = struct('X', 1*ones(3), 'Y', 2*ones(3), 'Z', 3*ones(3));          %
% end                                                                     %
% >> Y2 = cell(5,1);                                                      %
% >> for k = 1:5                                                          %
% Y2{k} = 2.5*ones(3);                                                    %
% end                                                                     %
% >> A = putfieldincell(A, 'Y2', Y2, 3);                                  %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% N/A                                                                     %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% putfieldincell
function [C] = putfieldincell(data, field, put, varargin)

% Determine the number of cells in the cell array 'data'.
n = length(data);

% Initialize the new cell array of structs.
C = data;

% Check if an index has been prescribed and that the index is valid.
if nargin == 4
    ind = varargin{1};
    nf  = numel(fieldnames(C{1}));
    if ind > 0 && floor(ind) == ind && ind <= nf+1
        flag = 1;
    end
end

if ~flag
    % If no index has been prescribed, place the new field at the end of
    % the field lists.
    for k = 1:n
        C{k}.(field) = put{k};
    end
else
    % If an index has been prescribed, place the new field at the
    % prescribed index and shift the remaining fields downward.
    for k = 1:n
        C{k}.(field) = put{k};
        if ind == 1
            perm = [ind, 1:nf];
        elseif ind == nf+1
            perm = 1:nf+1;
        else
            perm = [1:ind-1, nf+1, ind:nf];
        end
        C{k} = orderfields(C{k},perm);
    end
end

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
% * N/A.                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
