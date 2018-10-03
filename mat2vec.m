%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%           ARRANGE THE COLUMNS OF A MATRIX INTO A SINGLE COLUMN          %
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
% x = mat2vec(X);                                                         %
% x = mat2vec(X, 'nozeros');                                              %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function arranges the columns of a matrix into a single vector in  %
% sequential order.                                                       %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'x'          - 1D DOUBLE ARRAY                                          %
%              - Column vector to store the columns of a matrix 'X'.      %
% ----------------------------------------------------------------------- %
% 'X'          - 2D DOUBLE ARRAY                                          %
%              - Matrix to convert to a single column vector.             %
% ----------------------------------------------------------------------- %
% Options:                                                                %
% ----------------------------------------------------------------------- %
% 'nozeros'    - SPECIFIC STRING                                          %
%              - Option to ignore zero entries in the columnizing of 'X'  %
%                using 'nozeros'.                                         %
%              - Default: Not Active                                      %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Given a matrix A, write all the columns into a single vector.           %
%                                                                         %
% >> A = [1 4 7; 2 5 8; 3 6 9];                                           %
% >> a = mat2vec(A);                                                      %
%                                                                         %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% N/A                                                                     %
%                                                                         %
% Called in:                                                              %
% PP_CreateBlobs                                                          %
% PP_CreatePatches                                                        %
% PP_Write2TecPlot                                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% mat2vec
function [x] = mat2vec(X, varargin)

% Determine the size of the matrix 'X'.
[m,n] = size(X);

% Initialize the column vector to store all the columns of 'X'.
x = zeros(m*n,1);

if nargin == 2 && strcmpi(varargin{1}, 'nozeros')
    % Arrange the columns of 'X' sequentially in 'x'.
    c = 0;
    for j = 1:n
        for i = 1:m
            if X(i,j)
                c = c + 1;
                x(c) = X(i,j);
            else
                continue;
            end
        end
    end
    x(c+1:end) = [];
elseif nargin == 1
    % Arrange the columns of 'X' sequentially in 'x'.
    for k = 1:n
        x((k-1)*m+1:k*m) = X(:,k);
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
% Line(s) 63-68                                                           %
% * Note that the rows of a matrix could also be sequentially written in  %
%   a column vector by first taking its transpose, i.e. line 68 would be  %
%   written as:                                                           %
%   >> a = mat2vec(A.');                                                  %
% * Likewise, if a row vector output is required instead of a column      %
%   vector, one could simply take the transpose of the output.            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
