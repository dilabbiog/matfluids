%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            GEOMETRY TOOLBOX                             %
%                                                                         %
% idNodesub                                                               %
% Node Property Identification                                            %
% Examine geometry for subscripts of leading, interior and trailing nodes %
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
% [subl, subi, subt, subs] = idNodesub(geom);                             %
% [subl, subi, subt, subs] = idNodesub(geom, drc);                        %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Determine the subscripts of the leading, interior and trailing nodes of %
% an input array as well as any singular nodes. By default, subscripts    %
% are listed in the order they appear along the first dimension, since    %
% this is how linear indices are counted by MATLAB. Alternatively, the    %
% subscripts may be listed in the order they appear along a different     %
% direction.                                                              %
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
% 'geom'         LOGICAL/NUMERIC N-DIMENSIONAL ARRAY                      %
%              ~ Input array. Geometry or mask. Values are 1 within the   %
%                geometry and 0 outside the geometry.                     %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'drc'          POSITIVE INTEGER SCALAR                                  %
%                Default: 1                                               %
%              ~ Input scalar. Direction along which to list the          %
%                subscripts.                                              %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'subl'         CELL ARRAY                                               %
%              ~ Output cell array of size ndims(geom) by 1. Each cell k  %
%                contains the list of subscript k of leading nodes in the %
%                order they appear along a specified direction.           %
% ----------------------------------------------------------------------- %
% 'subi'         CELL ARRAY                                               %
%              ~ Output cell array of size ndims(geom) by 1. Each cell k  %
%                contains the list of subscript k of interior nodes in    %
%                the order they appear along a specified direction.       %
% ----------------------------------------------------------------------- %
% 'subt'         CELL ARRAY                                               %
%              ~ Output cell array of size ndims(geom) by 1. Each cell k  %
%                contains the list of subscript k of trailing nodes in    %
%                the order they appear along a specified direction.       %
% ----------------------------------------------------------------------- %
% 'subs'         CELL ARRAY                                               %
%              ~ Output cell array of size ndims(geom) by 1. Each cell k  %
%                contains the list of subscript k of singular nodes in    %
%                the order they appear along a specified direction.       %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Determine the leading, interior and trailing nodes for a geometry       %
% defined in one dimension.                                               %
%                                                                         %
% >> geom                     = [ones(20,1); zeros(10,1); ones(6,1)];     %
% >> [subl, subi, subt, subs] = idNodesub(geom);                          %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Determine the leading, interior and trailing nodes for a geometry       %
% defined in two dimensions.                                              %
%                                                                         %
% >> geom                     = zeros(51, 101);                           %
% >> geom(5:20,5:20)          = 1;                                        %
% >> geom(40:45,55:80)        = 1;                                        %
% >> [subl, subi, subt, subs] = idNodesub(geom, 1);                       %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 3                                                               %
%                                                                         %
% Determine the leading, interior and trailing nodes for a geometry       %
% defined in two dimensions. List the subscripts in the order they appear %
% along the second dimension.                                             %
%                                                                         %
% >> geom                     = zeros(51, 101);                           %
% >> geom(5:20,5:20)          = 1;                                        %
% >> geom(40:45,55:80)        = 1;                                        %
% >> [subl, subi, subt, subs] = idNodesub(geom, 2);                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [subl, subi, subt, subs] = idNodesub(geom, varargin)


%% PARSE INPUTS

% Input defaults.
default.drc = 1;

% Input checks.
check.geom = @(x) validateattributes(x,                                 ...
                  {'logical', 'numeric'},                               ...
                  {'nonempty', 'binary'});
check.drc  = @(x) validateattributes(x,                                 ...
                  {'logical', 'numeric'},                               ...
                  {'finite', 'positive', 'scalar'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'geom' , check.geom              );
addOptional ( hParser, 'drc'  , default.drc , check.drc );
parse(hParser, geom, varargin{:});
clear check default;

% Additional verifications.
nargoutchk(0,4);


%% INDEX SEARCH

% Determine the size and number of dimensions of the input array.
sz = size(geom);
nd = ndims(geom);

% Define the search direction and its length.
drc = hParser.Results.drc;
len = sz(drc);

% Permute the input array.
pm = 1:nd;
if drc > 1
    pm(drc) = [];
    pm      = [drc pm];
    geom    = permute(geom, pm);
end
szp = sz(pm);

% Strategy:
% 1) Find all nonzero elements of the input array in linear index form.
%    These represent all the nodes in the permuted geometry in the order
%    they appear moving down column by column.
     idxNZ = find(geom);
% 2) Compute the differences between successive nodes in "next - current"
%    form.
     id = [diff(idxNZ); 0];
% 3) All leading nodes and interior nodes are identified with +1. There is
%    also a false positive case identified with +1, namely, when a node is
%    at the end of the column (modulo len) and a subsequent node exists at
%    the beginning of the next column. This false positive is corrected in
%    step 4b.
% 4) Identify all possible trailing and singular nodes with 0. These
%    include:
%    a) All nodes identified with a value >1.
        id(id > 1) = 0;
%    b) All nodes bounded by the array.
        id(mod(idxNZ,len) == 0) = 0;
%    c) The last node. This node has already been set to 0 in step 2 to
%       have the differences in "next - current" form.
% 5) Compute the differences between the successive id values of the nodes,
%    this time in "current - previous" form. Note that the very first node
%    can either be a leading node (id = 1) or a singular node (id = 0). In
%    order to account for the first node in step 6a, when shifting the
%    differences down one index ("current - previous" form), use the actual
%    value of the first node.
     idd = [id(1); diff(id)];
% 6) All leading and trailing nodes can now be identified.
%    a) When a current node has id = 1 and a previous node has id = 0, the
%       current node must be a leading node and the difference will have a
%       value of +1. Distinguish these nodes using a value of 2.
        id(idd == 1) = 2;
%    b) When a current node has id = 0 and a previous node has id = 1, the
%       current node must be a trailing node and the difference will have a
%       value of -1. Distinguish these nodes using a value of 3.
        id(idd == -1) = 3;
%    c) When a current node has id = 0 and a previous node has id = 0, the
%       current node must be a singular node and the difference will have a
%       value of 0. These nodes are already identified by 0, id remains
%       unchanged.
%    d) When a current node has id = 1 and a previous node has id = 1, the
%       current node must be an interior node and the difference will have
%       a value of 0. These nodes are already identified by 1, id remains
%       unchanged.

% Define the leading, interior, trailing and singular nodes.
idxl = idxNZ(id == 2);
idxi = idxNZ(id == 1);
idxt = idxNZ(id == 3);
idxs = idxNZ(id == 0);
clear id idd idxNZ;

% Determine the leading node subscripts.
subl      = cell(nd,1);
[subl{:}] = ind2sub(szp, idxl);
clear idxl;

% Determine the interior node subscripts.
subi      = cell(nd,1);
[subi{:}] = ind2sub(szp, idxi);
clear idxi;

% Determine the trailing node subscripts.
subt      = cell(nd,1);
[subt{:}] = ind2sub(szp, idxt);
clear idxt;

% Determine the singular node subscripts.
subs      = cell(nd,1);
[subs{:}] = ind2sub(szp, idxs);
clear idxs;

% Determine the subscripts for the unpermuted input array.
if drc > 1
    subl = [subl(2:pm); subl(1); subl(pm+1:end)];
    subi = [subi(2:pm); subi(1); subi(pm+1:end)];
    subt = [subt(2:pm); subt(1); subt(pm+1:end)];
    subs = [subs(2:pm); subs(1); subs(pm+1:end)];
end
clear drc len nd pm sz szp;


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
% 2022/03/02 -- (GDL) Added missing support for singular nodes.           %
% 2022/02/28 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
