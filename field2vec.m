%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%          EXTRACT A STRUCT FIELD FROM A CELL ARRAY TO A 1D ARRAY         %
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
% C = field2vec(data, field);                                             %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function extracts a specified field with a single variable from a  %
% one-dimensional cell array of structs and outputs a column vector for   %
% that field.                                                             %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'C'          - 1D DOUBLE ARRAY                                          %
%              - One-dimensional array containing the extracted field     %
%                specified by the user.                                   %
% ----------------------------------------------------------------------- %
% 'data'       - 1D CELL                                                  %
%              - Cell array of structs with a field to be extracted.      %
% ----------------------------------------------------------------------- %
% 'field'      - STRING                                                   %
%              - The name of the field to be extracted.                   %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Given a cell array A of structs with fields 'X', 'Y', and 'Z', extract  %
% 'X' to have its own column vector.                                      %
%                                                                         %
% >> A = cell(5, 1);                                                      %
% >> for k = 1:5                                                          %
% A{k} = struct('X', k, 'Y', 2*k, 'Z', 3*k);                              %
% end                                                                     %
% >> X = field2vec(A, 'X');                                               %
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

%% field2vec
function [C] = field2vec(data, field)

% Determine the number of cells in the cell array 'data'.
n = length(data);

% Extract the specified field from the cell array of structs 'data'.
C = zeros(n, 1);
for k = 1:n
    C(k) = data{k}.(field);
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
% Line(s) 90                                                              %
% * Following MATLAB's recommendation: "Use dynamic fieldnames with       %
%   structures instead of GETFIELD", line 90 was selected over using      %
%   getfield(data{k}, field).                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
