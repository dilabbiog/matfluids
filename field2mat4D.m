%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%          EXTRACT A STRUCT FIELD FROM A CELL ARRAY TO A 4D ARRAY         %
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
% C = field2mat4D(data, field);                                           %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function extracts a specified field with a 3D matrix from a        %
% one-dimensional cell array of structs and outputs a 4D matrix for that  %
% field, with the fourth dimension representing the cell index.           %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'C'          - 4D DOUBLE ARRAY                                          %
%              - Four-dimensional array containing the extracted field    %
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
% 'X' to have its own 4D matrix.                                          %
%                                                                         %
% >> A = cell(5, 1);                                                      %
% >> for k = 1:5                                                          %
% A{k} = struct('X', k*ones(20,20,20), 'Y', 2*k*ones(20,20,20), 'Z', ...  %
%               3*k*ones(20,20,20));                                      %
% end                                                                     %
% >> X = field2mat4D(A, 'X');                                             %
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

%% field2mat4D
function [C] = field2mat4D(data, field)

% Determine the number of cells in the cell array 'data'.
n = length(data);

% Determine the number of rows, columns and 'pages'.
[r,c,p] = size(data{1}.(field));

% Extract the specified field from the cell array of structs 'data'.
C = zeros(r,c,p,n);
for k = 1:n
    C(:,:,:,k) = data{k}.(field);
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
% Line(s) 89,94                                                           %
% * Following MATLAB's recommendation: "Use dynamic fieldnames with       %
%   structures instead of GETFIELD", lines 89 and 94 were selected over   %
%   using getfield(data{k}, field).                                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
