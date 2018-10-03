%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                   TENSOR EIGENVALUES AND EIGENVECTORS                   %
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
% [E1, E2] = PP_TensorEig(TSR, opt);                                      %
% [E1, E2] = PP_TensorEig(TSR, opt, mask);                                %
% [E1, E2, E3] = PP_TensorEig(TSR, opt);                                  %
% [E1, E2, E3] = PP_TensorEig(TSR, opt, mask);                            %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function solve the tensor eigensystem associated with a two or     %
% three dimensional tensor with four or nine components respectively. The %
% function returns both the eigenvalues and eigenvectors of the tensor.   %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'E1'         - STRUCT                                                   %
%              - Struct containing the first eigenvalue in the 'eVal'     %
%                field and the components of the associated eigenvector   %
%                in the 'eVec1', 'eVec2' and (for 3D) 'eVec3' fields.     %
% ----------------------------------------------------------------------- %
% 'E2'         - STRUCT                                                   %
%              - Struct containing the second eigenvalue in the 'eVal'    %
%                field and the components of the associated eigenvector   %
%                in the 'eVec1', 'eVec2' and (for 3D) 'eVec3' fields.     %
% ----------------------------------------------------------------------- %
% 'E3'         - STRUCT                                                   %
%              - Struct containing the third eigenvalue in the 'eVal'     %
%                field and the components of the associated eigenvector   %
%                in the 'eVec1', 'eVec2' and (for 3D) 'eVec3' fields.     %
% ----------------------------------------------------------------------- %
% 'mask'       - 2D DOUBLE ARRAY                                          %
%              - Mask information for the tensor field. Values of 1       %
%                correspond to active nodes while values of 0 correspond  %
%                to blanked nodes.                                        %
% ----------------------------------------------------------------------- %
% 'opt'        - SPECIFIC STRING                                          %
%              - Option to premultiply the tensor by its transpose on the %
%                'left', on the 'right' or not at all using 'none'.       %
% ----------------------------------------------------------------------- %
% 'TSR'        - STRUCT                                                   %
%              - Two or three dimensional tensor with four or nine        %
%                components respectively. Each component represents a     %
%                field in a struct list. The order of the field list is   %
%                the order used to create a matrix elementwise by row for %
%                the eigensystem.                                         %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Solve the tensor eigensystem associated with the velocity gradient      %
% tensor of a time-dependent double gyre on the domain (x,y) = [0,2]x     %
% [0,1] with a constant grid spacing of 0.01 over the time interval       %
% [0,20] with time-step size 0.1. Use A = 0.1, epsilon = 0.25, and omega  %
% = 2*pi/10. Do not premultiply the tensor by its transpose.              %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> E1 = cell(length(t),1);                                              %
% >> E2 = cell(length(t),1);                                              %
% >> for k = 1:length(t)                                                  %
% [E1{k} E2{k}] = PP_TensorEig(VGT{k}, 'none');                           %
% end                                                                     %
% >> clear k;                                                             %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% N/A                                                                     %
%                                                                         %
% Called in:                                                              %
% PP_SwirlStrength                                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_TensorEig
function [varargout] = PP_TensorEig(TSR, opt, varargin)

fields = fieldnames(TSR);

if length(fields) == 4
    
    [ny,nx] = size(TSR.(fields{1}));
    
    varargout{1} = struct('eVal', zeros(ny,nx), 'eVec1', zeros(ny,nx), ...
                          'eVec2', zeros(ny,nx));
    varargout{2} = struct('eVal', zeros(ny,nx), 'eVec1', zeros(ny,nx), ...
                          'eVec2', zeros(ny,nx));
    
    switch nargin
        case 2
            if strcmpi(opt,'right')
                for i = 1:ny
                    for j = 1:nx
                        
                        C = [TSR.(fields{1})(i,j) TSR.(fields{2})(i,j);
                             TSR.(fields{3})(i,j) TSR.(fields{4})(i,j)];
                        C = C*C';
                        [eVecs, eVals] = eig(C);
                        [eVals, I] = sort(diag(eVals), 'ascend');
                        eVecs = eVecs(:,I.');
                        varargout{1}.eVal(i,j)  = eVals(1);
                        varargout{2}.eVal(i,j)  = eVals(2);
                        varargout{1}.eVec1(i,j) = eVecs(1,1);
                        varargout{1}.eVec2(i,j) = eVecs(2,1);
                        varargout{2}.eVec1(i,j) = eVecs(1,2);
                        varargout{2}.eVec2(i,j) = eVecs(2,2);
                        
                    end
                end
            elseif strcmpi(opt,'left')
                for i = 1:ny
                    for j = 1:nx
                        
                        C = [TSR.(fields{1})(i,j) TSR.(fields{2})(i,j);
                             TSR.(fields{3})(i,j) TSR.(fields{4})(i,j)];
                        C = C'*C;
                        [eVecs, eVals] = eig(C);
                        [eVals, I] = sort(diag(eVals), 'ascend');
                        eVecs = eVecs(:,I.');
                        varargout{1}.eVal(i,j)  = eVals(1);
                        varargout{2}.eVal(i,j)  = eVals(2);
                        varargout{1}.eVec1(i,j) = eVecs(1,1);
                        varargout{1}.eVec2(i,j) = eVecs(2,1);
                        varargout{2}.eVec1(i,j) = eVecs(1,2);
                        varargout{2}.eVec2(i,j) = eVecs(2,2);
                        
                    end
                end
            else
                for i = 1:ny
                    for j = 1:nx
                        
                        C = [TSR.(fields{1})(i,j) TSR.(fields{2})(i,j);
                             TSR.(fields{3})(i,j) TSR.(fields{4})(i,j)];
                        [eVecs, eVals] = eig(C);
                        [eVals, I] = sort(diag(eVals), 'ascend');
                        eVecs = eVecs(:,I.');
                        varargout{1}.eVal(i,j)  = eVals(1);
                        varargout{2}.eVal(i,j)  = eVals(2);
                        varargout{1}.eVec1(i,j) = eVecs(1,1);
                        varargout{1}.eVec2(i,j) = eVecs(2,1);
                        varargout{2}.eVec1(i,j) = eVecs(1,2);
                        varargout{2}.eVec2(i,j) = eVecs(2,2);
                        
                    end
                end
            end
        case 3
            if strcmpi(opt,'right')
                for i = 1:ny
                    for j = 1:nx
                        
                        if varargin{1}(i,j)
                            C = [TSR.(fields{1})(i,j) TSR.(fields{2})(i,j);
                                 TSR.(fields{3})(i,j) TSR.(fields{4})(i,j)];
                            C = C*C';
                            [eVecs, eVals] = eig(C);
                            [eVals, I] = sort(diag(eVals), 'ascend');
                            eVecs = eVecs(:,I.');
                            varargout{1}.eVal(i,j)  = eVals(1);
                            varargout{2}.eVal(i,j)  = eVals(2);
                            varargout{1}.eVec1(i,j) = eVecs(1,1);
                            varargout{1}.eVec2(i,j) = eVecs(2,1);
                            varargout{2}.eVec1(i,j) = eVecs(1,2);
                            varargout{2}.eVec2(i,j) = eVecs(2,2);
                        end
                        
                    end
                end
            elseif strcmpi(opt,'left')
                for i = 1:ny
                    for j = 1:nx
                        
                        if varargin{1}(i,j)
                            C = [TSR.(fields{1})(i,j) TSR.(fields{2})(i,j);
                                 TSR.(fields{3})(i,j) TSR.(fields{4})(i,j)];
                            C = C'*C;
                            [eVecs, eVals] = eig(C);
                            [eVals, I] = sort(diag(eVals), 'ascend');
                            eVecs = eVecs(:,I.');
                            varargout{1}.eVal(i,j)  = eVals(1);
                            varargout{2}.eVal(i,j)  = eVals(2);
                            varargout{1}.eVec1(i,j) = eVecs(1,1);
                            varargout{1}.eVec2(i,j) = eVecs(2,1);
                            varargout{2}.eVec1(i,j) = eVecs(1,2);
                            varargout{2}.eVec2(i,j) = eVecs(2,2);
                        end
                        
                    end
                end
            else
                for i = 1:ny
                    for j = 1:nx
                        
                        if varargin{1}(i,j)
                            C = [TSR.(fields{1})(i,j) TSR.(fields{2})(i,j);
                                 TSR.(fields{3})(i,j) TSR.(fields{4})(i,j)];
                            [eVecs, eVals] = eig(C);
                            [eVals, I] = sort(diag(eVals), 'ascend');
                            eVecs = eVecs(:,I.');
                            varargout{1}.eVal(i,j)  = eVals(1);
                            varargout{2}.eVal(i,j)  = eVals(2);
                            varargout{1}.eVec1(i,j) = eVecs(1,1);
                            varargout{1}.eVec2(i,j) = eVecs(2,1);
                            varargout{2}.eVec1(i,j) = eVecs(1,2);
                            varargout{2}.eVec2(i,j) = eVecs(2,2);
                        end
                        
                    end
                end
            end
    end
    
elseif length(fields) == 9
    
    [ny,nx,nz] = size(TSR.(fields{1}));
    
    varargout{1} = struct('eVal',  zeros(ny,nx,nz), ...
                          'eVec1', zeros(ny,nx,nz), ...
                          'eVec2', zeros(ny,nx,nz), ...
                          'eVec3', zeros(ny,nx,nz));
    varargout{2} = struct('eVal',  zeros(ny,nx,nz), ...
                          'eVec1', zeros(ny,nx,nz), ...
                          'eVec2', zeros(ny,nx,nz), ...
                          'eVec3', zeros(ny,nx,nz));
    varargout{3} = struct('eVal',  zeros(ny,nx,nz), ...
                          'eVec1', zeros(ny,nx,nz), ...
                          'eVec2', zeros(ny,nx,nz), ...
                          'eVec3', zeros(ny,nx,nz));
    
    switch nargin
        case 2
            if strcmpi(opt,'right')
                for k = 1:nz
                    for i = 1:ny
                        for j = 1:nx
                            
                            C = [TSR.(fields{1})(i,j,k) TSR.(fields{2})(i,j,k) TSR.(fields{3})(i,j,k);
                                 TSR.(fields{4})(i,j,k) TSR.(fields{5})(i,j,k) TSR.(fields{6})(i,j,k);
                                 TSR.(fields{7})(i,j,k) TSR.(fields{8})(i,j,k) TSR.(fields{9})(i,j,k)];
                            C = C*C';
                            [eVecs, eVals] = eig(C);
                            [eVals, I] = sort(diag(eVals), 'ascend');
                            eVecs = eVecs(:,I.');
                            varargout{1}.eVal(i,j,k)  = eVals(1);
                            varargout{2}.eVal(i,j,k)  = eVals(2);
                            varargout{3}.eVal(i,j,k)  = eVals(2);
                            varargout{1}.eVec1(i,j,k) = eVecs(1,1);
                            varargout{1}.eVec2(i,j,k) = eVecs(2,1);
                            varargout{1}.eVec3(i,j,k) = eVecs(3,1);
                            varargout{2}.eVec1(i,j,k) = eVecs(1,2);
                            varargout{2}.eVec2(i,j,k) = eVecs(2,2);
                            varargout{2}.eVec3(i,j,k) = eVecs(3,2);
                            varargout{3}.eVec1(i,j,k) = eVecs(1,3);
                            varargout{3}.eVec2(i,j,k) = eVecs(2,3);
                            varargout{3}.eVec3(i,j,k) = eVecs(3,3);
                            
                        end
                    end
                end
            elseif strcmpi(opt,'left')
                for k = 1:nz
                    for i = 1:ny
                        for j = 1:nx
                            
                            C = [TSR.(fields{1})(i,j,k) TSR.(fields{2})(i,j,k) TSR.(fields{3})(i,j,k);
                                 TSR.(fields{4})(i,j,k) TSR.(fields{5})(i,j,k) TSR.(fields{6})(i,j,k);
                                 TSR.(fields{7})(i,j,k) TSR.(fields{8})(i,j,k) TSR.(fields{9})(i,j,k)];
                            C = C'*C;
                            [eVecs, eVals] = eig(C);
                            [eVals, I] = sort(diag(eVals), 'ascend');
                            eVecs = eVecs(:,I.');
                            varargout{1}.eVal(i,j,k)  = eVals(1);
                            varargout{2}.eVal(i,j,k)  = eVals(2);
                            varargout{3}.eVal(i,j,k)  = eVals(2);
                            varargout{1}.eVec1(i,j,k) = eVecs(1,1);
                            varargout{1}.eVec2(i,j,k) = eVecs(2,1);
                            varargout{1}.eVec3(i,j,k) = eVecs(3,1);
                            varargout{2}.eVec1(i,j,k) = eVecs(1,2);
                            varargout{2}.eVec2(i,j,k) = eVecs(2,2);
                            varargout{2}.eVec3(i,j,k) = eVecs(3,2);
                            varargout{3}.eVec1(i,j,k) = eVecs(1,3);
                            varargout{3}.eVec2(i,j,k) = eVecs(2,3);
                            varargout{3}.eVec3(i,j,k) = eVecs(3,3);
                            
                        end
                    end
                end
            else
                for k = 1:nz
                    for i = 1:ny
                        for j = 1:nx
                            
                            C = [TSR.(fields{1})(i,j,k) TSR.(fields{2})(i,j,k) TSR.(fields{3})(i,j,k);
                                 TSR.(fields{4})(i,j,k) TSR.(fields{5})(i,j,k) TSR.(fields{6})(i,j,k);
                                 TSR.(fields{7})(i,j,k) TSR.(fields{8})(i,j,k) TSR.(fields{9})(i,j,k)];
                            [eVecs, eVals] = eig(C);
                            [eVals, I] = sort(diag(eVals), 'ascend');
                            eVecs = eVecs(:,I.');
                            varargout{1}.eVal(i,j,k)  = eVals(1);
                            varargout{2}.eVal(i,j,k)  = eVals(2);
                            varargout{3}.eVal(i,j,k)  = eVals(2);
                            varargout{1}.eVec1(i,j,k) = eVecs(1,1);
                            varargout{1}.eVec2(i,j,k) = eVecs(2,1);
                            varargout{1}.eVec3(i,j,k) = eVecs(3,1);
                            varargout{2}.eVec1(i,j,k) = eVecs(1,2);
                            varargout{2}.eVec2(i,j,k) = eVecs(2,2);
                            varargout{2}.eVec3(i,j,k) = eVecs(3,2);
                            varargout{3}.eVec1(i,j,k) = eVecs(1,3);
                            varargout{3}.eVec2(i,j,k) = eVecs(2,3);
                            varargout{3}.eVec3(i,j,k) = eVecs(3,3);
                            
                        end
                    end
                end
            end
        case 3
            if strcmpi(opt,'right')
                for k = 1:nz
                    for i = 1:ny
                        for j = 1:nx
                            
                            if varargin{1}(i,j,k)
                                C = [TSR.(fields{1})(i,j,k) TSR.(fields{2})(i,j,k) TSR.(fields{3})(i,j,k);
                                     TSR.(fields{4})(i,j,k) TSR.(fields{5})(i,j,k) TSR.(fields{6})(i,j,k);
                                     TSR.(fields{7})(i,j,k) TSR.(fields{8})(i,j,k) TSR.(fields{9})(i,j,k)];
                                C = C*C';
                                [eVecs, eVals] = eig(C);
                                [eVals, I] = sort(diag(eVals), 'ascend');
                                eVecs = eVecs(:,I.');
                                varargout{1}.eVal(i,j,k)  = eVals(1);
                                varargout{2}.eVal(i,j,k)  = eVals(2);
                                varargout{3}.eVal(i,j,k)  = eVals(2);
                                varargout{1}.eVec1(i,j,k) = eVecs(1,1);
                                varargout{1}.eVec2(i,j,k) = eVecs(2,1);
                                varargout{1}.eVec3(i,j,k) = eVecs(3,1);
                                varargout{2}.eVec1(i,j,k) = eVecs(1,2);
                                varargout{2}.eVec2(i,j,k) = eVecs(2,2);
                                varargout{2}.eVec3(i,j,k) = eVecs(3,2);
                                varargout{3}.eVec1(i,j,k) = eVecs(1,3);
                                varargout{3}.eVec2(i,j,k) = eVecs(2,3);
                                varargout{3}.eVec3(i,j,k) = eVecs(3,3);
                            end
                            
                        end
                    end
                end
            elseif strcmpi(opt,'left')
                for k = 1:nz
                    for i = 1:ny
                        for j = 1:nx
                            
                            if varargin{1}(i,j,k)
                                C = [TSR.(fields{1})(i,j,k) TSR.(fields{2})(i,j,k) TSR.(fields{3})(i,j,k);
                                     TSR.(fields{4})(i,j,k) TSR.(fields{5})(i,j,k) TSR.(fields{6})(i,j,k);
                                     TSR.(fields{7})(i,j,k) TSR.(fields{8})(i,j,k) TSR.(fields{9})(i,j,k)];
                                C = C'*C;
                                [eVecs, eVals] = eig(C);
                                [eVals, I] = sort(diag(eVals), 'ascend');
                                eVecs = eVecs(:,I.');
                                varargout{1}.eVal(i,j,k)  = eVals(1);
                                varargout{2}.eVal(i,j,k)  = eVals(2);
                                varargout{3}.eVal(i,j,k)  = eVals(2);
                                varargout{1}.eVec1(i,j,k) = eVecs(1,1);
                                varargout{1}.eVec2(i,j,k) = eVecs(2,1);
                                varargout{1}.eVec3(i,j,k) = eVecs(3,1);
                                varargout{2}.eVec1(i,j,k) = eVecs(1,2);
                                varargout{2}.eVec2(i,j,k) = eVecs(2,2);
                                varargout{2}.eVec3(i,j,k) = eVecs(3,2);
                                varargout{3}.eVec1(i,j,k) = eVecs(1,3);
                                varargout{3}.eVec2(i,j,k) = eVecs(2,3);
                                varargout{3}.eVec3(i,j,k) = eVecs(3,3);
                            end
                            
                        end
                    end
                end
            else
                for k = 1:nz
                    for i = 1:ny
                        for j = 1:nx
                            
                            if varargin{1}(i,j,k)
                                C = [TSR.(fields{1})(i,j,k) TSR.(fields{2})(i,j,k) TSR.(fields{3})(i,j,k);
                                     TSR.(fields{4})(i,j,k) TSR.(fields{5})(i,j,k) TSR.(fields{6})(i,j,k);
                                     TSR.(fields{7})(i,j,k) TSR.(fields{8})(i,j,k) TSR.(fields{9})(i,j,k)];
                                [eVecs, eVals] = eig(C);
                                [eVals, I] = sort(diag(eVals), 'ascend');
                                eVecs = eVecs(:,I.');
                                varargout{1}.eVal(i,j,k)  = eVals(1);
                                varargout{2}.eVal(i,j,k)  = eVals(2);
                                varargout{3}.eVal(i,j,k)  = eVals(2);
                                varargout{1}.eVec1(i,j,k) = eVecs(1,1);
                                varargout{1}.eVec2(i,j,k) = eVecs(2,1);
                                varargout{1}.eVec3(i,j,k) = eVecs(3,1);
                                varargout{2}.eVec1(i,j,k) = eVecs(1,2);
                                varargout{2}.eVec2(i,j,k) = eVecs(2,2);
                                varargout{2}.eVec3(i,j,k) = eVecs(3,2);
                                varargout{3}.eVec1(i,j,k) = eVecs(1,3);
                                varargout{3}.eVec2(i,j,k) = eVecs(2,3);
                                varargout{3}.eVec3(i,j,k) = eVecs(3,3);
                            end
                            
                        end
                    end
                end
            end
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
