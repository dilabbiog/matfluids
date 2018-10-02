%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%             CONVERT A CELL ARRAY OF MATRICES INTO A 3D ARRAY            %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Giuseppe Di Labbio                                                      %
% Department of Mechanical, Industrial & Aerospace Engineering            %
% Concordia University Montréal, Canada                                   %
%                                                                         %
% Last Update: October 2nd, 2018 by Giuseppe Di Labbio                    %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% SYNTAX                                                                  %
%                                                                         %
% C = cell2mat3D(data);                                                   %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function converts a one-dimensional cell array of matrices to a 3D %
% matrix with the third dimension representing the cell index.            %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'C'          - 3D DOUBLE ARRAY                                          %
%              - Three-dimensional array containing the matrices of a     %
%                cell array along the third dimension. The index of the   %
%                third dimension is equivalent to the cell index.         %
% ----------------------------------------------------------------------- %
% 'data'       - 1D CELL, ELEMENTS: 2D DOUBLE ARRAYS                      %
%              - One-dimensional cell array of matrices.                  %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Transform a cell array of matrices A into a 3D matrix.                  %
%                                                                         %
% >> A = cell(5, 1);                                                      %
% >> for k = 1:5, A{k} = k*ones(5); end                                   %
% >> X = cell2mat3D(A);                                                   %
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

%% cell2mat3D
function [C] = cell2mat3D(data)

%% Parse the inputs.
if iscell(data)
    if ndims(data) == 2 && min(size(data)) == 1
        for k = 1:length(data)
            if isnumeric(data{k}) || islogical(data{k})
                if ~ismatrix(squeeze(data{k}))
                    error('InputError:isMatrix', ...
                         ['One or more elements of the cell array has ' ...
                          'more than two dimensions (after squeezing).']);
                end
                
                if k > 1
                    tmp = size(squeeze(data{k})) ...
                        - size(squeeze(data{k-1}));
                    if max(tmp)
                        error('InputError:sizeMismatch', ...
                             ['The sizes of two or more matrices in ' ...
                              'the cell array do not match.']);
                    end
                end
                
            else
                error('InputError:isNumeric', ...
                     ['One or more elements of the cell array is not ' ...
                      'numeric.']);
            end
        end
    else
        error('InputError:is1Dcell', ...
              'The cell array has more than one dimension.');
    end
else
    error('InputError:isCellArray', ...
          'Input is not a one-dimensional cell array.');
end
clear tmp;

%% Determine the number of cells in the cell array 'data'.
n = length(data);

%% Determine the number of rows and columns.
[r,c] = size(squeeze(data{1}));

%% Convert the cell array into a three-dimensional array.
C = zeros(r,c,n);
for k = 1:n
    C(:,:,k) = squeeze(data{k});
end
clear c k n r;

%% %%%%%%%%%%%%%%%%%%%%%%%%% SUPPRESS MESSAGES %%%%%%%%%%%%%%%%%%%%%%%%% %%

%#ok<*ISMAT>
% Line(s) 63
% Message(s)
% * When checking if a variable is a matrix consider using ISMATRIX.
% Reason(s)
% * The variable being inspected is a cell array, not a matrix.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Line(s) N/A                                                             %
% * N/A.                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
