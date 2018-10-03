%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            FINITE DIFFERENCES                           %
%              FOURTH ORDER HYBRID COMPACT RICHARDSON SCHEME              %
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
% DRV = FD_CR4(VEC);                                                      %
% DRV = FD_CR4(VEC, 'mask');                                              %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the spatial derivatives of velocity for 2D or 3D %
% Cartesian data. The derivatives are computed using the noise-optimized, %
% fourth-order hybrid compact Richardson scheme of Ali Etebari & Pavlos   %
% P. Vlachos [1] (which makes use of the schemes in [2]). The function    %
% switches to lower-order schemes in the cases where insufficent points   %
% are present.                                                            %
%                                                                         %
% References:                                                             %
% [1] Etebari, A. & Vlachos, P. P. (2005). Improvements on the accuracy   %
%     of derivative estimation from DPIV velocity measurements.           %
%     Experiments in Fluids, 39, 1040-1050.                               %
% [2] Lele, S. K. (1992). Compact finite difference schemes with          %
%     spectral-like resolution. Journal of Computational Physics, 103,    %
%     16-42.                                                              %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'DRV'        - STRUCT OF 2D OR 3D DOUBLE ARRAYS                         %
%              - Derivatives of the velocity components U, V and W in X,  %
%                Y and Z.                                                 %
% ----------------------------------------------------------------------- %
% 'VEC'        - STRUCT OF 2D OR 3D DOUBLE ARRAYS                         %
%              - Velocity field data containing the mask 'C', Cartesian   %
%                coordinates 'X', 'Y' and 'Z', and velocity field         %
%                components 'U', 'V' and 'W'.                             %
% ----------------------------------------------------------------------- %
% Options:                                                                %
% ----------------------------------------------------------------------- %
% 'mask'       - SPECIFIC STRING                                          %
%              - Option to use a pre-defined mask in VEC.C. If a mask is  %
%                available, the calculation will be considerably sped up  %
%                and the boundary nodes will be treated correctly.        %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Compute the derivatives of the double gyre flow with constants A = 0.1, %
% epsilon = 0.25 and omega = 2*pi/10 at t = 3. Use a domain of (x,y) =    %
% [0,2]x[0,1]. Compare the derivatives with their exact values.           %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t    = 3;                                                            %
% >> A    = 0.1;                                                          %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC,DRV] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                  %
% >> VGT = FD_CR4(VEC{1});                                                %
% >> Eux = max(max(DRV{1}.UX - VGT.UX));                                  %
% >> Euy = max(max(DRV{1}.UY - VGT.UY));                                  %
% >> Evx = max(max(DRV{1}.VX - VGT.VX));                                  %
% >> Evy = max(max(DRV{1}.VY - VGT.VY));                                  %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% FD_COMP2                                                                %
% FD_COMP3                                                                %
% FD_COMP4                                                                %
% FD_EXO2                                                                 %
% FD_SCH                                                                  %
% TDMA                                                                    %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FD_CR4
function [DRV] = FD_CR4(VEC, varargin)

if nargin == 1
    if length(fieldnames(VEC)) == 5
        
        % Determine the number of data points in x (nx) and y (ny).
        [ny, nx] = size(VEC.U);
        
        % Compute the grid spacing in the x (dx) and y (dy) directions.
        dx = VEC.X(1,2) - VEC.X(1,1);
        dy = VEC.Y(2,1) - VEC.Y(1,1);
        
        % Define the noise-optimized constants in the Richardson scheme.
        A  = 1239;
        A1 = 272;
        A2 = 1036;
        A4 = -69;
        
        % Declare a struct holding the derivatives of the functions U and V
        % in X and Y.
        DRV = struct('UX', 0, 'UY', 0, 'VX', 0, 'VY', 0);
        
        % Initialize the reduced-grid x derivative matrices.
        ux1 = zeros(ny,nx);
        ux2 = zeros(ny,nx);
        ux4 = zeros(ny,nx);
        vx1 = zeros(ny,nx);
        vx2 = zeros(ny,nx);
        vx4 = zeros(ny,nx);
        
        for i = 1:ny
            
            % Compute du/dx for the full grid (dx).
            ux1(i,:)      = FD_SCH(VEC.U(i,:), dx, 'COMP4').';
            
            % Compute du/dx for alternative grids skipping 1 node (2*dx).
            ux2(i,1:2:nx) = FD_SCH(VEC.U(i,1:2:nx), 2*dx, 'COMP4').';
            ux2(i,2:2:nx) = FD_SCH(VEC.U(i,2:2:nx), 2*dx, 'COMP4').';
            
            % Compute du/dx for alternative grids skipping 3 nodes (4*dx).
            ux4(i,1:4:nx) = FD_SCH(VEC.U(i,1:4:nx), 4*dx, 'COMP4').';
            ux4(i,2:4:nx) = FD_SCH(VEC.U(i,2:4:nx), 4*dx, 'COMP4').';
            ux4(i,3:4:nx) = FD_SCH(VEC.U(i,3:4:nx), 4*dx, 'COMP4').';
            ux4(i,4:4:nx) = FD_SCH(VEC.U(i,4:4:nx), 4*dx, 'COMP4').';
            
            % Compute dv/dx for the full grid (dx).
            vx1(i,:)      = FD_SCH(VEC.V(i,:), dx, 'COMP4').';
            
            % Compute dv/dx for alternative grids skipping 1 node (2*dx).
            vx2(i,1:2:nx) = FD_SCH(VEC.V(i,1:2:nx), 2*dx, 'COMP4').';
            vx2(i,2:2:nx) = FD_SCH(VEC.V(i,2:2:nx), 2*dx, 'COMP4').';
            
            % Compute dv/dx for alternative grids skipping 3 nodes (4*dx).
            vx4(i,1:4:nx) = FD_SCH(VEC.V(i,1:4:nx), 4*dx, 'COMP4').';
            vx4(i,2:4:nx) = FD_SCH(VEC.V(i,2:4:nx), 4*dx, 'COMP4').';
            vx4(i,3:4:nx) = FD_SCH(VEC.V(i,3:4:nx), 4*dx, 'COMP4').';
            vx4(i,4:4:nx) = FD_SCH(VEC.V(i,4:4:nx), 4*dx, 'COMP4').';
            
        end
        
        % Compute du/dx and dv/dx using the noise-optimized Richardson
        % scheme.
        DRV.UX = (1/A)*(A1*ux1 + A2*ux2 + A4*ux4);
        DRV.VX = (1/A)*(A1*vx1 + A2*vx2 + A4*vx4);
        
        % Initialize the reduced-grid y derivative matrices.
        uy1 = zeros(ny,nx);
        uy2 = zeros(ny,nx);
        uy4 = zeros(ny,nx);
        vy1 = zeros(ny,nx);
        vy2 = zeros(ny,nx);
        vy4 = zeros(ny,nx);
        
        for j = 1:nx
            
            % Compute du/dy for the full grid (dy).
            uy1(:,j)      = FD_SCH(VEC.U(:,j), dy, 'COMP4');
            
            % Compute du/dy for alternative grids skipping 1 node (2*dy).
            uy2(1:2:ny,j) = FD_SCH(VEC.U(1:2:ny,j), 2*dy, 'COMP4');
            uy2(2:2:ny,j) = FD_SCH(VEC.U(2:2:ny,j), 2*dy, 'COMP4');
            
            % Compute du/dy for alternative grids skipping 3 nodes (4*dy).
            uy4(1:4:ny,j) = FD_SCH(VEC.U(1:4:ny,j), 4*dy, 'COMP4');
            uy4(2:4:ny,j) = FD_SCH(VEC.U(2:4:ny,j), 4*dy, 'COMP4');
            uy4(3:4:ny,j) = FD_SCH(VEC.U(3:4:ny,j), 4*dy, 'COMP4');
            uy4(4:4:ny,j) = FD_SCH(VEC.U(4:4:ny,j), 4*dy, 'COMP4');
            
            % Compute dv/dy for the full grid (dy).
            vy1(:,j)      = FD_SCH(VEC.V(:,j), dy, 'COMP4');
            
            % Compute dv/dy for alternative grids skipping 1 node (2*dy).
            vy2(1:2:ny,j) = FD_SCH(VEC.V(1:2:ny,j), 2*dy, 'COMP4');
            vy2(2:2:ny,j) = FD_SCH(VEC.V(2:2:ny,j), 2*dy, 'COMP4');
            
            % Compute dv/dy for alternative grids skipping 3 nodes (4*dy).
            vy4(1:4:ny,j) = FD_SCH(VEC.V(1:4:ny,j), 4*dy, 'COMP4');
            vy4(2:4:ny,j) = FD_SCH(VEC.V(2:4:ny,j), 4*dy, 'COMP4');
            vy4(3:4:ny,j) = FD_SCH(VEC.V(3:4:ny,j), 4*dy, 'COMP4');
            vy4(4:4:ny,j) = FD_SCH(VEC.V(4:4:ny,j), 4*dy, 'COMP4');
            
        end
        
        % Compute du/dy and dv/dy using the noise-optimized Richardson
        % scheme.
        DRV.UY = (1/A)*(A1*uy1 + A2*uy2 + A4*uy4);
        DRV.VY = (1/A)*(A1*vy1 + A2*vy2 + A4*vy4);
        
    elseif length(fieldnames(VEC)) == 6
        
        % Determine the number of data points in x (nx) and y (ny).
        [ny, nx] = size(VEC.U);
        
        % Compute the grid spacing in the x (dx) and y (dy) directions.
        dx = VEC.X(1,2) - VEC.X(1,1);
        dy = VEC.Y(2,1) - VEC.Y(1,1);
        
        % Define the noise-optimized constants in the Richardson scheme.
        A  = 1239;
        A1 = 272;
        A2 = 1036;
        A4 = -69;
        
        % Declare a struct holding the derivatives of the functions U, V
        % and W in X and Y.
        DRV = struct('UX', zeros(ny,nx), 'UY', zeros(ny,nx), ...
                     'VX', zeros(ny,nx), 'VY', zeros(ny,nx), ...
                     'WX', zeros(ny,nx), 'WY', zeros(ny,nx));
        
        % Initialize the reduced-grid x derivative matrices.
        ux1 = zeros(ny,nx);
        ux2 = zeros(ny,nx);
        ux4 = zeros(ny,nx);
        vx1 = zeros(ny,nx);
        vx2 = zeros(ny,nx);
        vx4 = zeros(ny,nx);
        wx1 = zeros(ny,nx);
        wx2 = zeros(ny,nx);
        wx4 = zeros(ny,nx);
        
        for i = 1:ny
            
            % Compute du/dx for the full grid (dx).
            ux1(i,:)      = FD_SCH(VEC.U(i,:), dx, 'COMP4').';
            
            % Compute du/dx for alternative grids skipping 1 node
            % (2*dx).
            ux2(i,1:2:nx) = FD_SCH(VEC.U(i,1:2:nx), 2*dx, 'COMP4').';
            ux2(i,2:2:nx) = FD_SCH(VEC.U(i,2:2:nx), 2*dx, 'COMP4').';
            
            % Compute du/dx for alternative grids skipping 3 nodes
            % (4*dx).
            ux4(i,1:4:nx) = FD_SCH(VEC.U(i,1:4:nx), 4*dx, 'COMP4').';
            ux4(i,2:4:nx) = FD_SCH(VEC.U(i,2:4:nx), 4*dx, 'COMP4').';
            ux4(i,3:4:nx) = FD_SCH(VEC.U(i,3:4:nx), 4*dx, 'COMP4').';
            ux4(i,4:4:nx) = FD_SCH(VEC.U(i,4:4:nx), 4*dx, 'COMP4').';
            
            % Compute dv/dx for the full grid (dx).
            vx1(i,:)      = FD_SCH(VEC.V(i,:), dx, 'COMP4').';
            
            % Compute dv/dx for alternative grids skipping 1 node
            % (2*dx).
            vx2(i,1:2:nx) = FD_SCH(VEC.V(i,1:2:nx), 2*dx, 'COMP4').';
            vx2(i,2:2:nx) = FD_SCH(VEC.V(i,2:2:nx), 2*dx, 'COMP4').';
            
            % Compute dv/dx for alternative grids skipping 3 nodes
            % (4*dx).
            vx4(i,1:4:nx) = FD_SCH(VEC.V(i,1:4:nx), 4*dx, 'COMP4').';
            vx4(i,2:4:nx) = FD_SCH(VEC.V(i,2:4:nx), 4*dx, 'COMP4').';
            vx4(i,3:4:nx) = FD_SCH(VEC.V(i,3:4:nx), 4*dx, 'COMP4').';
            vx4(i,4:4:nx) = FD_SCH(VEC.V(i,4:4:nx), 4*dx, 'COMP4').';
            
            % Compute dw/dx for the full grid (dx).
            wx1(i,:)      = FD_SCH(VEC.W(i,:), dx, 'COMP4').';
            
            % Compute dw/dx for alternative grids skipping 1 node
            % (2*dx).
            wx2(i,1:2:nx) = FD_SCH(VEC.W(i,1:2:nx), 2*dx, 'COMP4').';
            wx2(i,2:2:nx) = FD_SCH(VEC.W(i,2:2:nx), 2*dx, 'COMP4').';
            
            % Compute dw/dx for alternative grids skipping 3 nodes
            % (4*dx).
            wx4(i,1:4:nx) = FD_SCH(VEC.W(i,1:4:nx), 4*dx, 'COMP4').';
            wx4(i,2:4:nx) = FD_SCH(VEC.W(i,2:4:nx), 4*dx, 'COMP4').';
            wx4(i,3:4:nx) = FD_SCH(VEC.W(i,3:4:nx), 4*dx, 'COMP4').';
            wx4(i,4:4:nx) = FD_SCH(VEC.W(i,4:4:nx), 4*dx, 'COMP4').';
            
        end
        
        % Compute du/dx, dv/dx and dw/dx using the noise-optimized
        % Richardson scheme.
        DRV.UX = (1/A)*(A1*ux1 + A2*ux2 + A4*ux4);
        DRV.VX = (1/A)*(A1*vx1 + A2*vx2 + A4*vx4);
        DRV.WX = (1/A)*(A1*wx1 + A2*wx2 + A4*wx4);
        
        % Initialize the reduced-grid y derivative matrices.
        uy1 = zeros(ny,nx);
        uy2 = zeros(ny,nx);
        uy4 = zeros(ny,nx);
        vy1 = zeros(ny,nx);
        vy2 = zeros(ny,nx);
        vy4 = zeros(ny,nx);
        wy1 = zeros(ny,nx);
        wy2 = zeros(ny,nx);
        wy4 = zeros(ny,nx);
        
        for j = 1:nx
            
            % Compute du/dy for the full grid (dy).
            uy1(:,j)      = FD_SCH(VEC.U(:,j), dy, 'COMP4');
            
            % Compute du/dy for alternative grids skipping 1 node
            % (2*dy).
            uy2(1:2:ny,j) = FD_SCH(VEC.U(1:2:ny,j), 2*dy, 'COMP4');
            uy2(2:2:ny,j) = FD_SCH(VEC.U(2:2:ny,j), 2*dy, 'COMP4');
            
            % Compute du/dy for alternative grids skipping 3 nodes
            % (4*dy).
            uy4(1:4:ny,j) = FD_SCH(VEC.U(1:4:ny,j), 4*dy, 'COMP4');
            uy4(2:4:ny,j) = FD_SCH(VEC.U(2:4:ny,j), 4*dy, 'COMP4');
            uy4(3:4:ny,j) = FD_SCH(VEC.U(3:4:ny,j), 4*dy, 'COMP4');
            uy4(4:4:ny,j) = FD_SCH(VEC.U(4:4:ny,j), 4*dy, 'COMP4');
            
            % Compute dv/dy for the full grid (dy).
            vy1(:,j)      = FD_SCH(VEC.V(:,j), dy, 'COMP4');
            
            % Compute dv/dy for alternative grids skipping 1 node
            % (2*dy).
            vy2(1:2:ny,j) = FD_SCH(VEC.V(1:2:ny,j), 2*dy, 'COMP4');
            vy2(2:2:ny,j) = FD_SCH(VEC.V(2:2:ny,j), 2*dy, 'COMP4');
            
            % Compute dv/dy for alternative grids skipping 3 nodes
            % (4*dy).
            vy4(1:4:ny,j) = FD_SCH(VEC.V(1:4:ny,j), 4*dy, 'COMP4');
            vy4(2:4:ny,j) = FD_SCH(VEC.V(2:4:ny,j), 4*dy, 'COMP4');
            vy4(3:4:ny,j) = FD_SCH(VEC.V(3:4:ny,j), 4*dy, 'COMP4');
            vy4(4:4:ny,j) = FD_SCH(VEC.V(4:4:ny,j), 4*dy, 'COMP4');
            
            % Compute dw/dy for the full grid (dy).
            wy1(:,j)      = FD_SCH(VEC.W(:,j), dy, 'COMP4');
            
            % Compute dw/dy for alternative grids skipping 1 node
            % (2*dy).
            wy2(1:2:ny,j) = FD_SCH(VEC.W(1:2:ny,j), 2*dy, 'COMP4');
            wy2(2:2:ny,j) = FD_SCH(VEC.W(2:2:ny,j), 2*dy, 'COMP4');
            
            % Compute dw/dy for alternative grids skipping 3 nodes
            % (4*dy).
            wy4(1:4:ny,j) = FD_SCH(VEC.W(1:4:ny,j), 4*dy, 'COMP4');
            wy4(2:4:ny,j) = FD_SCH(VEC.W(2:4:ny,j), 4*dy, 'COMP4');
            wy4(3:4:ny,j) = FD_SCH(VEC.W(3:4:ny,j), 4*dy, 'COMP4');
            wy4(4:4:ny,j) = FD_SCH(VEC.W(4:4:ny,j), 4*dy, 'COMP4');
            
        end
        
        % Compute du/dy, dv/dy and dw/dy using the noise-optimized
        % Richardson scheme.
        DRV.UY = (1/A)*(A1*uy1 + A2*uy2 + A4*uy4);
        DRV.VY = (1/A)*(A1*vy1 + A2*vy2 + A4*vy4);
        DRV.WY = (1/A)*(A1*wy1 + A2*wy2 + A4*wy4);
            
    elseif length(fieldnames(VEC)) == 7
        
        % Determine the number of data points in x (nx), y (ny) and z (nz).
        [ny, nx, nz] = size(VEC.U);
        
        % Compute the grid spacing in the x (dx), y (dy) and z (dz)
        % directions.
        dx = VEC.X(1,2) - VEC.X(1,1);
        dy = VEC.Y(2,1) - VEC.Y(1,1);
        dz = VEC.Z(2)   - VEC.Z(1);
        
        % Define the noise-optimized constants in the Richardson scheme.
        A  = 1239;
        A1 = 272;
        A2 = 1036;
        A4 = -69;
        
        % Declare a struct holding the derivatives of the functions U, V
        % and W in X, Y and Z.
        DRV = struct('UX', zeros(ny,nx,nz), 'UY', zeros(ny,nx,nz), ...
                     'UZ', zeros(ny,nx,nz), 'VX', zeros(ny,nx,nz), ...
                     'VY', zeros(ny,nx,nz), 'VZ', zeros(ny,nx,nz), ...
                     'WX', zeros(ny,nx,nz), 'WY', zeros(ny,nx,nz), ...
                     'WZ', zeros(ny,nx,nz));
        
        % Initialize the reduced-grid x derivative matrices.
        ux1 = zeros(ny,nx);
        ux2 = zeros(ny,nx);
        ux4 = zeros(ny,nx);
        vx1 = zeros(ny,nx);
        vx2 = zeros(ny,nx);
        vx4 = zeros(ny,nx);
        wx1 = zeros(ny,nx);
        wx2 = zeros(ny,nx);
        wx4 = zeros(ny,nx);
        
        for k = 1:nz
            for i = 1:ny
                
                % Compute du/dx for the full grid (dx).
                ux1(i,:)      = FD_SCH(VEC.U(i,:,k), dx, 'COMP4').';
                
                % Compute du/dx for alternative grids skipping 1 node
                % (2*dx).
                ux2(i,1:2:nx) = FD_SCH(VEC.U(i,1:2:nx,k), 2*dx, 'COMP4').';
                ux2(i,2:2:nx) = FD_SCH(VEC.U(i,2:2:nx,k), 2*dx, 'COMP4').';
                
                % Compute du/dx for alternative grids skipping 3 nodes
                % (4*dx).
                ux4(i,1:4:nx) = FD_SCH(VEC.U(i,1:4:nx,k), 4*dx, 'COMP4').';
                ux4(i,2:4:nx) = FD_SCH(VEC.U(i,2:4:nx,k), 4*dx, 'COMP4').';
                ux4(i,3:4:nx) = FD_SCH(VEC.U(i,3:4:nx,k), 4*dx, 'COMP4').';
                ux4(i,4:4:nx) = FD_SCH(VEC.U(i,4:4:nx,k), 4*dx, 'COMP4').';
                
                % Compute dv/dx for the full grid (dx).
                vx1(i,:)      = FD_SCH(VEC.V(i,:,k), dx, 'COMP4').';
                
                % Compute dv/dx for alternative grids skipping 1 node
                % (2*dx).
                vx2(i,1:2:nx) = FD_SCH(VEC.V(i,1:2:nx,k), 2*dx, 'COMP4').';
                vx2(i,2:2:nx) = FD_SCH(VEC.V(i,2:2:nx,k), 2*dx, 'COMP4').';
                
                % Compute dv/dx for alternative grids skipping 3 nodes
                % (4*dx).
                vx4(i,1:4:nx) = FD_SCH(VEC.V(i,1:4:nx,k), 4*dx, 'COMP4').';
                vx4(i,2:4:nx) = FD_SCH(VEC.V(i,2:4:nx,k), 4*dx, 'COMP4').';
                vx4(i,3:4:nx) = FD_SCH(VEC.V(i,3:4:nx,k), 4*dx, 'COMP4').';
                vx4(i,4:4:nx) = FD_SCH(VEC.V(i,4:4:nx,k), 4*dx, 'COMP4').';
                
                % Compute dw/dx for the full grid (dx).
                wx1(i,:)      = FD_SCH(VEC.W(i,:,k), dx, 'COMP4').';
                
                % Compute dw/dx for alternative grids skipping 1 node
                % (2*dx).
                wx2(i,1:2:nx) = FD_SCH(VEC.W(i,1:2:nx,k), 2*dx, 'COMP4').';
                wx2(i,2:2:nx) = FD_SCH(VEC.W(i,2:2:nx,k), 2*dx, 'COMP4').';
                
                % Compute dw/dx for alternative grids skipping 3 nodes
                % (4*dx).
                wx4(i,1:4:nx) = FD_SCH(VEC.W(i,1:4:nx,k), 4*dx, 'COMP4').';
                wx4(i,2:4:nx) = FD_SCH(VEC.W(i,2:4:nx,k), 4*dx, 'COMP4').';
                wx4(i,3:4:nx) = FD_SCH(VEC.W(i,3:4:nx,k), 4*dx, 'COMP4').';
                wx4(i,4:4:nx) = FD_SCH(VEC.W(i,4:4:nx,k), 4*dx, 'COMP4').';
                
            end
            
            % Compute du/dx, dv/dx and dw/dx using the noise-optimized
            % Richardson scheme.
            DRV.UX(:,:,k) = (1/A)*(A1*ux1 + A2*ux2 + A4*ux4);
            DRV.VX(:,:,k) = (1/A)*(A1*vx1 + A2*vx2 + A4*vx4);
            DRV.WX(:,:,k) = (1/A)*(A1*wx1 + A2*wx2 + A4*wx4);
            
        end
        
        % Initialize the reduced-grid y derivative matrices.
        uy1 = zeros(ny,nx);
        uy2 = zeros(ny,nx);
        uy4 = zeros(ny,nx);
        vy1 = zeros(ny,nx);
        vy2 = zeros(ny,nx);
        vy4 = zeros(ny,nx);
        wy1 = zeros(ny,nx);
        wy2 = zeros(ny,nx);
        wy4 = zeros(ny,nx);
        
        for k = 1:nz
            for j = 1:nx
                
                % Compute du/dy for the full grid (dy).
                uy1(:,j)      = FD_SCH(VEC.U(:,j,k), dy, 'COMP4');
                
                % Compute du/dy for alternative grids skipping 1 node
                % (2*dy).
                uy2(1:2:ny,j) = FD_SCH(VEC.U(1:2:ny,j,k), 2*dy, 'COMP4');
                uy2(2:2:ny,j) = FD_SCH(VEC.U(2:2:ny,j,k), 2*dy, 'COMP4');
                
                % Compute du/dy for alternative grids skipping 3 nodes
                % (4*dy).
                uy4(1:4:ny,j) = FD_SCH(VEC.U(1:4:ny,j,k), 4*dy, 'COMP4');
                uy4(2:4:ny,j) = FD_SCH(VEC.U(2:4:ny,j,k), 4*dy, 'COMP4');
                uy4(3:4:ny,j) = FD_SCH(VEC.U(3:4:ny,j,k), 4*dy, 'COMP4');
                uy4(4:4:ny,j) = FD_SCH(VEC.U(4:4:ny,j,k), 4*dy, 'COMP4');
                
                % Compute dv/dy for the full grid (dy).
                vy1(:,j)      = FD_SCH(VEC.V(:,j,k), dy, 'COMP4');
                
                % Compute dv/dy for alternative grids skipping 1 node
                % (2*dy).
                vy2(1:2:ny,j) = FD_SCH(VEC.V(1:2:ny,j,k), 2*dy, 'COMP4');
                vy2(2:2:ny,j) = FD_SCH(VEC.V(2:2:ny,j,k), 2*dy, 'COMP4');
                
                % Compute dv/dy for alternative grids skipping 3 nodes
                % (4*dy).
                vy4(1:4:ny,j) = FD_SCH(VEC.V(1:4:ny,j,k), 4*dy, 'COMP4');
                vy4(2:4:ny,j) = FD_SCH(VEC.V(2:4:ny,j,k), 4*dy, 'COMP4');
                vy4(3:4:ny,j) = FD_SCH(VEC.V(3:4:ny,j,k), 4*dy, 'COMP4');
                vy4(4:4:ny,j) = FD_SCH(VEC.V(4:4:ny,j,k), 4*dy, 'COMP4');
                
                % Compute dw/dy for the full grid (dy).
                wy1(:,j)      = FD_SCH(VEC.W(:,j,k), dy, 'COMP4');
                
                % Compute dw/dy for alternative grids skipping 1 node
                % (2*dy).
                wy2(1:2:ny,j) = FD_SCH(VEC.W(1:2:ny,j,k), 2*dy, 'COMP4');
                wy2(2:2:ny,j) = FD_SCH(VEC.W(2:2:ny,j,k), 2*dy, 'COMP4');
                
                % Compute dw/dy for alternative grids skipping 3 nodes
                % (4*dy).
                wy4(1:4:ny,j) = FD_SCH(VEC.W(1:4:ny,j,k), 4*dy, 'COMP4');
                wy4(2:4:ny,j) = FD_SCH(VEC.W(2:4:ny,j,k), 4*dy, 'COMP4');
                wy4(3:4:ny,j) = FD_SCH(VEC.W(3:4:ny,j,k), 4*dy, 'COMP4');
                wy4(4:4:ny,j) = FD_SCH(VEC.W(4:4:ny,j,k), 4*dy, 'COMP4');
                
            end
            
            % Compute du/dy, dv/dy and dw/dy using the noise-optimized
            % Richardson scheme.
            DRV.UY(:,:,k) = (1/A)*(A1*uy1 + A2*uy2 + A4*uy4);
            DRV.VY(:,:,k) = (1/A)*(A1*vy1 + A2*vy2 + A4*vy4);
            DRV.WY(:,:,k) = (1/A)*(A1*wy1 + A2*wy2 + A4*wy4);
            
        end
        
        % Initialize the reduced-grid z derivative matrices.
        uz1 = zeros(ny,nx,nz);
        uz2 = zeros(ny,nx,nz);
        uz4 = zeros(ny,nx,nz);
        vz1 = zeros(ny,nx,nz);
        vz2 = zeros(ny,nx,nz);
        vz4 = zeros(ny,nx,nz);
        wz1 = zeros(ny,nx,nz);
        wz2 = zeros(ny,nx,nz);
        wz4 = zeros(ny,nx,nz);
        
        for i = 1:ny
            for j = 1:nx
                
                % Compute du/dz for the full grid (dz).
                uz1(i,j,:)      = FD_SCH(VEC.U(i,j,:), dz, 'COMP4');
                
                % Compute du/dz for alternative grids skipping 1 node
                % (2*dz).
                uz2(i,j,1:2:nz) = FD_SCH(VEC.U(i,j,1:2:nz), 2*dz, 'COMP4');
                uz2(i,j,2:2:nz) = FD_SCH(VEC.U(i,j,2:2:nz), 2*dz, 'COMP4');
                
                % Compute du/dz for alternative grids skipping 3 nodes
                % (4*dz).
                uz4(i,j,1:4:nz) = FD_SCH(VEC.U(i,j,1:4:nz), 4*dz, 'COMP4');
                uz4(i,j,2:4:nz) = FD_SCH(VEC.U(i,j,2:4:nz), 4*dz, 'COMP4');
                uz4(i,j,3:4:nz) = FD_SCH(VEC.U(i,j,3:4:nz), 4*dz, 'COMP4');
                uz4(i,j,4:4:nz) = FD_SCH(VEC.U(i,j,4:4:nz), 4*dz, 'COMP4');
                
                % Compute dv/dz for the full grid (dz).
                vz1(i,j,:)      = FD_SCH(VEC.V(i,j,:), dz, 'COMP4');
                
                % Compute dv/dz for alternative grids skipping 1 node
                % (2*dz).
                vz2(i,j,1:2:nz) = FD_SCH(VEC.V(i,j,1:2:nz), 2*dz, 'COMP4');
                vz2(i,j,2:2:nz) = FD_SCH(VEC.V(i,j,2:2:nz), 2*dz, 'COMP4');
                
                % Compute dv/dz for alternative grids skipping 3 nodes
                % (4*dz).
                vz4(i,j,1:4:nz) = FD_SCH(VEC.V(i,j,1:4:nz), 4*dz, 'COMP4');
                vz4(i,j,2:4:nz) = FD_SCH(VEC.V(i,j,2:4:nz), 4*dz, 'COMP4');
                vz4(i,j,3:4:nz) = FD_SCH(VEC.V(i,j,3:4:nz), 4*dz, 'COMP4');
                vz4(i,j,4:4:nz) = FD_SCH(VEC.V(i,j,4:4:nz), 4*dz, 'COMP4');
                
                % Compute dw/dz for the full grid (dz).
                wz1(i,j,:)      = FD_SCH(VEC.W(i,j,:), dz, 'COMP4');
                
                % Compute dw/dz for alternative grids skipping 1 node
                % (2*dz).
                wz2(i,j,1:2:nz) = FD_SCH(VEC.W(i,j,1:2:nz), 2*dz, 'COMP4');
                wz2(i,j,2:2:nz) = FD_SCH(VEC.W(i,j,2:2:nz), 2*dz, 'COMP4');
                
                % Compute dw/dz for alternative grids skipping 3 nodes
                % (4*dz).
                wz4(i,j,1:4:nz) = FD_SCH(VEC.W(i,j,1:4:nz), 4*dz, 'COMP4');
                wz4(i,j,2:4:nz) = FD_SCH(VEC.W(i,j,2:4:nz), 4*dz, 'COMP4');
                wz4(i,j,3:4:nz) = FD_SCH(VEC.W(i,j,3:4:nz), 4*dz, 'COMP4');
                wz4(i,j,4:4:nz) = FD_SCH(VEC.W(i,j,4:4:nz), 4*dz, 'COMP4');
                
            end
        end
        
        % Compute du/dz, dv/dz and dw/dz using the noise-optimized
        % Richardson scheme.
        DRV.UZ = (1/A)*(A1*uz1 + A2*uz2 + A4*uz4);
        DRV.VZ = (1/A)*(A1*vz1 + A2*vz2 + A4*vz4);
        DRV.WZ = (1/A)*(A1*wz1 + A2*wz2 + A4*wz4);
    end
        
elseif nargin == 2 && strcmpi(varargin{1}, 'mask')
    if length(fieldnames(VEC)) == 5
        
        % Determine the number of data points in x (nx) and y (ny).
        [ny, nx] = size(VEC.U);
        
        % Compute the grid spacing in the x (dx) and y (dy) directions.
        dx = VEC.X(1,2) - VEC.X(1,1);
        dy = VEC.Y(2,1) - VEC.Y(1,1);
        
        % Define the noise-optimized constants in the Richardson scheme.
        A  = 1239;
        A1 = 272;
        A2 = 1036;
        A4 = -69;
        
        % Declare a struct holding the derivatives of the functions U and V
        % in X and Y.
        DRV = struct('UX', 0, 'UY', 0, 'VX', 0, 'VY', 0);
        
        % Initialize the reduced-grid x derivative matrices.
        ux1 = zeros(ny,nx);
        ux2 = zeros(ny,nx);
        ux4 = zeros(ny,nx);
        vx1 = zeros(ny,nx);
        vx2 = zeros(ny,nx);
        vx4 = zeros(ny,nx);
        
        for i = 1:ny
            
            % Compute du/dx for the full grid (dx).
            ux1(i,:)      = FD_SCH(VEC.U(i,:), dx, 'COMP4', VEC.C(i,:)).';
            
            % Compute du/dx for alternative grids skipping 1 node (2*dx).
            ux2(i,1:2:nx) = FD_SCH(VEC.U(i,1:2:nx), 2*dx, ...
                                   'COMP4', VEC.C(i,1:2:nx)).';
            ux2(i,2:2:nx) = FD_SCH(VEC.U(i,2:2:nx), 2*dx, ...
                                   'COMP4', VEC.C(i,2:2:nx)).';
            
            % Compute du/dx for alternative grids skipping 3 nodes (4*dx).
            ux4(i,1:4:nx) = FD_SCH(VEC.U(i,1:4:nx), 4*dx, ...
                                   'COMP4', VEC.C(i,1:4:nx)).';
            ux4(i,2:4:nx) = FD_SCH(VEC.U(i,2:4:nx), 4*dx, ...
                                   'COMP4', VEC.C(i,2:4:nx)).';
            ux4(i,3:4:nx) = FD_SCH(VEC.U(i,3:4:nx), 4*dx, ...
                                   'COMP4', VEC.C(i,3:4:nx)).';
            ux4(i,4:4:nx) = FD_SCH(VEC.U(i,4:4:nx), 4*dx, ...
                                   'COMP4', VEC.C(i,4:4:nx)).';
            
            % Compute dv/dx for the full grid (dx).
            vx1(i,:)      = FD_SCH(VEC.V(i,:), dx, 'COMP4', VEC.C(i,:)).';
            
            % Compute dv/dx for alternative grids skipping 1 node (2*dx).
            vx2(i,1:2:nx) = FD_SCH(VEC.V(i,1:2:nx), 2*dx, ...
                                   'COMP4', VEC.C(i,1:2:nx)).';
            vx2(i,2:2:nx) = FD_SCH(VEC.V(i,2:2:nx), 2*dx, ...
                                   'COMP4', VEC.C(i,2:2:nx)).';
            
            % Compute dv/dx for alternative grids skipping 3 nodes (4*dx).
            vx4(i,1:4:nx) = FD_SCH(VEC.V(i,1:4:nx), 4*dx, ...
                                   'COMP4', VEC.C(i,1:4:nx)).';
            vx4(i,2:4:nx) = FD_SCH(VEC.V(i,2:4:nx), 4*dx, ...
                                   'COMP4', VEC.C(i,2:4:nx)).';
            vx4(i,3:4:nx) = FD_SCH(VEC.V(i,3:4:nx), 4*dx, ...
                                   'COMP4', VEC.C(i,3:4:nx)).';
            vx4(i,4:4:nx) = FD_SCH(VEC.V(i,4:4:nx), 4*dx, ...
                                   'COMP4', VEC.C(i,4:4:nx)).';
            
        end
        
        % Compute du/dx and dv/dx using the noise-optimized Richardson
        % scheme.
        DRV.UX = (1/A)*(A1*ux1 + A2*ux2 + A4*ux4);
        DRV.VX = (1/A)*(A1*vx1 + A2*vx2 + A4*vx4);
        
        % Initialize the reduced-grid y derivative matrices.
        uy1 = zeros(ny,nx);
        uy2 = zeros(ny,nx);
        uy4 = zeros(ny,nx);
        vy1 = zeros(ny,nx);
        vy2 = zeros(ny,nx);
        vy4 = zeros(ny,nx);
        
        for j = 1:nx
            
            % Compute du/dy for the full grid (dy).
            uy1(:,j)      = FD_SCH(VEC.U(:,j), dy, 'COMP4', VEC.C(:,j));
            
            % Compute du/dy for alternative grids skipping 1 node (2*dy).
            uy2(1:2:ny,j) = FD_SCH(VEC.U(1:2:ny,j), 2*dy, ...
                                   'COMP4', VEC.C(1:2:ny,j));
            uy2(2:2:ny,j) = FD_SCH(VEC.U(2:2:ny,j), 2*dy, ...
                                   'COMP4', VEC.C(2:2:ny,j));
            
            % Compute du/dy for alternative grids skipping 3 nodes (4*dy).
            uy4(1:4:ny,j) = FD_SCH(VEC.U(1:4:ny,j), 4*dy, ...
                                   'COMP4', VEC.C(1:4:ny,j));
            uy4(2:4:ny,j) = FD_SCH(VEC.U(2:4:ny,j), 4*dy, ...
                                   'COMP4', VEC.C(2:4:ny,j));
            uy4(3:4:ny,j) = FD_SCH(VEC.U(3:4:ny,j), 4*dy, ...
                                   'COMP4', VEC.C(3:4:ny,j));
            uy4(4:4:ny,j) = FD_SCH(VEC.U(4:4:ny,j), 4*dy, ...
                                   'COMP4', VEC.C(4:4:ny,j));
            
            % Compute dv/dy for the full grid (dy).
            vy1(:,j)      = FD_SCH(VEC.V(:,j), dy, 'COMP4', VEC.C(:,j));
            
            % Compute dv/dy for alternative grids skipping 1 node (2*dy).
            vy2(1:2:ny,j) = FD_SCH(VEC.V(1:2:ny,j), 2*dy, ...
                                   'COMP4', VEC.C(1:2:ny,j));
            vy2(2:2:ny,j) = FD_SCH(VEC.V(2:2:ny,j), 2*dy, ...
                                   'COMP4', VEC.C(2:2:ny,j));
            
            % Compute dv/dy for alternative grids skipping 3 nodes (4*dy).
            vy4(1:4:ny,j) = FD_SCH(VEC.V(1:4:ny,j), 4*dy, ...
                                   'COMP4', VEC.C(1:4:ny,j));
            vy4(2:4:ny,j) = FD_SCH(VEC.V(2:4:ny,j), 4*dy, ...
                                   'COMP4', VEC.C(2:4:ny,j));
            vy4(3:4:ny,j) = FD_SCH(VEC.V(3:4:ny,j), 4*dy, ...
                                   'COMP4', VEC.C(3:4:ny,j));
            vy4(4:4:ny,j) = FD_SCH(VEC.V(4:4:ny,j), 4*dy, ...
                                   'COMP4', VEC.C(4:4:ny,j));
            
        end
        
        % Compute du/dy and dv/dy using the noise-optimized Richardson
        % scheme.
        DRV.UY = (1/A)*(A1*uy1 + A2*uy2 + A4*uy4);
        DRV.VY = (1/A)*(A1*vy1 + A2*vy2 + A4*vy4);
    
    elseif length(fieldnames(VEC)) == 6
        
        % Determine the number of data points in x (nx) and y (ny).
        [ny, nx] = size(VEC.U);
        
        % Compute the grid spacing in the x (dx) and y (dy) directions.
        dx = VEC.X(1,2) - VEC.X(1,1);
        dy = VEC.Y(2,1) - VEC.Y(1,1);
        
        % Define the noise-optimized constants in the Richardson scheme.
        A  = 1239;
        A1 = 272;
        A2 = 1036;
        A4 = -69;
        
        % Declare a struct holding the derivatives of the functions U, V
        % and W in X and Y.
        DRV = struct('UX', zeros(ny,nx), 'UY', zeros(ny,nx), ...
                     'VX', zeros(ny,nx), 'VY', zeros(ny,nx), ...
                     'WX', zeros(ny,nx), 'WY', zeros(ny,nx));
        
        % Initialize the reduced-grid x derivative matrices.
        ux1 = zeros(ny,nx);
        ux2 = zeros(ny,nx);
        ux4 = zeros(ny,nx);
        vx1 = zeros(ny,nx);
        vx2 = zeros(ny,nx);
        vx4 = zeros(ny,nx);
        wx1 = zeros(ny,nx);
        wx2 = zeros(ny,nx);
        wx4 = zeros(ny,nx);
        
        for i = 1:ny
            
            % Compute du/dx for the full grid (dx).
            ux1(i,:)      = FD_SCH(VEC.U(i,:), dx, 'COMP4', VEC.C(i,:)).';
            
            % Compute du/dx for alternative grids skipping 1 node
            % (2*dx).
            ux2(i,1:2:nx) = FD_SCH(VEC.U(i,1:2:nx), 2*dx, 'COMP4', ...
                                   VEC.C(i,1:2:nx)).';
            ux2(i,2:2:nx) = FD_SCH(VEC.U(i,2:2:nx), 2*dx, 'COMP4', ...
                                   VEC.C(i,2:2:nx)).';
            
            % Compute du/dx for alternative grids skipping 3 nodes
            % (4*dx).
            ux4(i,1:4:nx) = FD_SCH(VEC.U(i,1:4:nx), 4*dx, 'COMP4', ...
                                   VEC.C(i,1:4:nx)).';
            ux4(i,2:4:nx) = FD_SCH(VEC.U(i,2:4:nx), 4*dx, 'COMP4', ...
                                   VEC.C(i,2:4:nx)).';
            ux4(i,3:4:nx) = FD_SCH(VEC.U(i,3:4:nx), 4*dx, 'COMP4', ...
                                   VEC.C(i,3:4:nx)).';
            ux4(i,4:4:nx) = FD_SCH(VEC.U(i,4:4:nx), 4*dx, 'COMP4', ...
                                   VEC.C(i,4:4:nx)).';
            
            % Compute dv/dx for the full grid (dx).
            vx1(i,:)      = FD_SCH(VEC.V(i,:), dx, 'COMP4', VEC.C(i,:)).';
            
            % Compute dv/dx for alternative grids skipping 1 node
            % (2*dx).
            vx2(i,1:2:nx) = FD_SCH(VEC.V(i,1:2:nx), 2*dx, 'COMP4', ...
                                   VEC.C(i,1:2:nx)).';
            vx2(i,2:2:nx) = FD_SCH(VEC.V(i,2:2:nx), 2*dx, 'COMP4', ...
                                   VEC.C(i,2:2:nx)).';
            
            % Compute dv/dx for alternative grids skipping 3 nodes
            % (4*dx).
            vx4(i,1:4:nx) = FD_SCH(VEC.V(i,1:4:nx), 4*dx, 'COMP4', ...
                                   VEC.C(i,1:4:nx)).';
            vx4(i,2:4:nx) = FD_SCH(VEC.V(i,2:4:nx), 4*dx, 'COMP4', ...
                                   VEC.C(i,2:4:nx)).';
            vx4(i,3:4:nx) = FD_SCH(VEC.V(i,3:4:nx), 4*dx, 'COMP4', ...
                                   VEC.C(i,3:4:nx)).';
            vx4(i,4:4:nx) = FD_SCH(VEC.V(i,4:4:nx), 4*dx, 'COMP4', ...
                                   VEC.C(i,4:4:nx)).';
            
            % Compute dw/dx for the full grid (dx).
            wx1(i,:)      = FD_SCH(VEC.W(i,:), dx, 'COMP4', VEC.C(i,:)).';
            
            % Compute dw/dx for alternative grids skipping 1 node
            % (2*dx).
            wx2(i,1:2:nx) = FD_SCH(VEC.W(i,1:2:nx), 2*dx, 'COMP4', ...
                                   VEC.C(i,1:2:nx)).';
            wx2(i,2:2:nx) = FD_SCH(VEC.W(i,2:2:nx), 2*dx, 'COMP4', ...
                                   VEC.C(i,2:2:nx)).';
            
            % Compute dw/dx for alternative grids skipping 3 nodes
            % (4*dx).
            wx4(i,1:4:nx) = FD_SCH(VEC.W(i,1:4:nx), 4*dx, 'COMP4', ...
                                   VEC.C(i,1:4:nx)).';
            wx4(i,2:4:nx) = FD_SCH(VEC.W(i,2:4:nx), 4*dx, 'COMP4', ...
                                   VEC.C(i,2:4:nx)).';
            wx4(i,3:4:nx) = FD_SCH(VEC.W(i,3:4:nx), 4*dx, 'COMP4', ...
                                   VEC.C(i,3:4:nx)).';
            wx4(i,4:4:nx) = FD_SCH(VEC.W(i,4:4:nx), 4*dx, 'COMP4', ...
                                   VEC.C(i,4:4:nx)).';
            
        end
        
        % Compute du/dx, dv/dx and dw/dx using the noise-optimized
        % Richardson scheme.
        DRV.UX = (1/A)*(A1*ux1 + A2*ux2 + A4*ux4);
        DRV.VX = (1/A)*(A1*vx1 + A2*vx2 + A4*vx4);
        DRV.WX = (1/A)*(A1*wx1 + A2*wx2 + A4*wx4);
        
        % Initialize the reduced-grid y derivative matrices.
        uy1 = zeros(ny,nx);
        uy2 = zeros(ny,nx);
        uy4 = zeros(ny,nx);
        vy1 = zeros(ny,nx);
        vy2 = zeros(ny,nx);
        vy4 = zeros(ny,nx);
        wy1 = zeros(ny,nx);
        wy2 = zeros(ny,nx);
        wy4 = zeros(ny,nx);
        
        for j = 1:nx
            
            % Compute du/dy for the full grid (dy).
            uy1(:,j)      = FD_SCH(VEC.U(:,j), dy, 'COMP4', VEC.C(:,j));
            
            % Compute du/dy for alternative grids skipping 1 node
            % (2*dy).
            uy2(1:2:ny,j) = FD_SCH(VEC.U(1:2:ny,j), 2*dy, 'COMP4', ...
                                   VEC.C(1:2:ny,j));
            uy2(2:2:ny,j) = FD_SCH(VEC.U(2:2:ny,j), 2*dy, 'COMP4', ...
                                   VEC.C(2:2:ny,j));
            
            % Compute du/dy for alternative grids skipping 3 nodes
            % (4*dy).
            uy4(1:4:ny,j) = FD_SCH(VEC.U(1:4:ny,j), 4*dy, 'COMP4', ...
                                   VEC.C(1:4:ny,j));
            uy4(2:4:ny,j) = FD_SCH(VEC.U(2:4:ny,j), 4*dy, 'COMP4', ...
                                   VEC.C(2:4:ny,j));
            uy4(3:4:ny,j) = FD_SCH(VEC.U(3:4:ny,j), 4*dy, 'COMP4', ...
                                   VEC.C(3:4:ny,j));
            uy4(4:4:ny,j) = FD_SCH(VEC.U(4:4:ny,j), 4*dy, 'COMP4', ...
                                   VEC.C(4:4:ny,j));
            
            % Compute dv/dy for the full grid (dy).
            vy1(:,j)      = FD_SCH(VEC.V(:,j), dy, 'COMP4', VEC.C(:,j));
            
            % Compute dv/dy for alternative grids skipping 1 node
            % (2*dy).
            vy2(1:2:ny,j) = FD_SCH(VEC.V(1:2:ny,j), 2*dy, 'COMP4', ...
                                   VEC.C(1:2:ny,j));
            vy2(2:2:ny,j) = FD_SCH(VEC.V(2:2:ny,j), 2*dy, 'COMP4', ...
                                   VEC.C(2:2:ny,j));
            
            % Compute dv/dy for alternative grids skipping 3 nodes
            % (4*dy).
            vy4(1:4:ny,j) = FD_SCH(VEC.V(1:4:ny,j), 4*dy, 'COMP4', ...
                                   VEC.C(1:4:ny,j));
            vy4(2:4:ny,j) = FD_SCH(VEC.V(2:4:ny,j), 4*dy, 'COMP4', ...
                                   VEC.C(2:4:ny,j));
            vy4(3:4:ny,j) = FD_SCH(VEC.V(3:4:ny,j), 4*dy, 'COMP4', ...
                                   VEC.C(3:4:ny,j));
            vy4(4:4:ny,j) = FD_SCH(VEC.V(4:4:ny,j), 4*dy, 'COMP4', ...
                                   VEC.C(4:4:ny,j));
            
            % Compute dw/dy for the full grid (dy).
            wy1(:,j)      = FD_SCH(VEC.W(:,j), dy, 'COMP4', VEC.C(:,j));
            
            % Compute dw/dy for alternative grids skipping 1 node
            % (2*dy).
            wy2(1:2:ny,j) = FD_SCH(VEC.W(1:2:ny,j), 2*dy, 'COMP4', ...
                                   VEC.C(1:2:ny,j));
            wy2(2:2:ny,j) = FD_SCH(VEC.W(2:2:ny,j), 2*dy, 'COMP4', ...
                                   VEC.C(2:2:ny,j));
            
            % Compute dw/dy for alternative grids skipping 3 nodes
            % (4*dy).
            wy4(1:4:ny,j) = FD_SCH(VEC.W(1:4:ny,j), 4*dy, 'COMP4', ...
                                   VEC.C(1:4:ny,j));
            wy4(2:4:ny,j) = FD_SCH(VEC.W(2:4:ny,j), 4*dy, 'COMP4', ...
                                   VEC.C(2:4:ny,j));
            wy4(3:4:ny,j) = FD_SCH(VEC.W(3:4:ny,j), 4*dy, 'COMP4', ...
                                   VEC.C(3:4:ny,j));
            wy4(4:4:ny,j) = FD_SCH(VEC.W(4:4:ny,j), 4*dy, 'COMP4', ...
                                   VEC.C(4:4:ny,j));
            
        end
        
        % Compute du/dy, dv/dy and dw/dy using the noise-optimized
        % Richardson scheme.
        DRV.UY = (1/A)*(A1*uy1 + A2*uy2 + A4*uy4);
        DRV.VY = (1/A)*(A1*vy1 + A2*vy2 + A4*vy4);
        DRV.WY = (1/A)*(A1*wy1 + A2*wy2 + A4*wy4);
        
    elseif length(fieldnames(VEC)) == 7
        
        % Determine the number of data points in x (nx), y (ny) and z (nz).
        [ny, nx, nz] = size(VEC.U);
        
        % Compute the grid spacing in the x (dx), y (dy) and z (dz)
        % directions.
        dx = VEC.X(1,2) - VEC.X(1,1);
        dy = VEC.Y(2,1) - VEC.Y(1,1);
        dz = VEC.Z(2)   - VEC.Z(1);
        
        % Define the noise-optimized constants in the Richardson scheme.
        A  = 1239;
        A1 = 272;
        A2 = 1036;
        A4 = -69;
        
        % Declare a struct holding the derivatives of the functions U, V
        % and W in X, Y and Z.
        DRV = struct('UX', zeros(ny,nx,nz), 'UY', zeros(ny,nx,nz), ...
                     'UZ', zeros(ny,nx,nz), 'VX', zeros(ny,nx,nz), ...
                     'VY', zeros(ny,nx,nz), 'VZ', zeros(ny,nx,nz), ...
                     'WX', zeros(ny,nx,nz), 'WY', zeros(ny,nx,nz), ...
                     'WZ', zeros(ny,nx,nz));
        
        % Initialize the reduced-grid x derivative matrices.
        ux1 = zeros(ny,nx);
        ux2 = zeros(ny,nx);
        ux4 = zeros(ny,nx);
        vx1 = zeros(ny,nx);
        vx2 = zeros(ny,nx);
        vx4 = zeros(ny,nx);
        wx1 = zeros(ny,nx);
        wx2 = zeros(ny,nx);
        wx4 = zeros(ny,nx);
        
        for k = 1:nz
            for i = 1:ny
                
                % Compute du/dx for the full grid (dx).
                ux1(i,:)      = FD_SCH(VEC.U(i,:,k), dx, 'COMP4', ...
                                       VEC.C(i,:,k)).';
                
                % Compute du/dx for alternative grids skipping 1 node
                % (2*dx).
                ux2(i,1:2:nx) = FD_SCH(VEC.U(i,1:2:nx,k), 2*dx, ...
                                       'COMP4', VEC.C(i,1:2:nx,k)).';
                ux2(i,2:2:nx) = FD_SCH(VEC.U(i,2:2:nx,k), 2*dx, ...
                                       'COMP4', VEC.C(i,2:2:nx,k)).';
                
                % Compute du/dx for alternative grids skipping 3 nodes
                % (4*dx).
                ux4(i,1:4:nx) = FD_SCH(VEC.U(i,1:4:nx,k), 4*dx, ...
                                       'COMP4', VEC.C(i,1:4:nx,k)).';
                ux4(i,2:4:nx) = FD_SCH(VEC.U(i,2:4:nx,k), 4*dx, ...
                                       'COMP4', VEC.C(i,2:4:nx,k)).';
                ux4(i,3:4:nx) = FD_SCH(VEC.U(i,3:4:nx,k), 4*dx, ...
                                       'COMP4', VEC.C(i,3:4:nx,k)).';
                ux4(i,4:4:nx) = FD_SCH(VEC.U(i,4:4:nx,k), 4*dx, ...
                                       'COMP4', VEC.C(i,4:4:nx,k)).';
                
                % Compute dv/dx for the full grid (dx).
                vx1(i,:)      = FD_SCH(VEC.V(i,:,k), dx, 'COMP4', ...
                                       VEC.C(i,:,k)).';
                
                % Compute dv/dx for alternative grids skipping 1 node
                % (2*dx).
                vx2(i,1:2:nx) = FD_SCH(VEC.V(i,1:2:nx,k), 2*dx, ...
                                       'COMP4', VEC.C(i,1:2:nx)).';
                vx2(i,2:2:nx) = FD_SCH(VEC.V(i,2:2:nx,k), 2*dx, ...
                                       'COMP4', VEC.C(i,2:2:nx)).';
                
                % Compute dv/dx for alternative grids skipping 3 nodes
                % (4*dx).
                vx4(i,1:4:nx) = FD_SCH(VEC.V(i,1:4:nx,k), 4*dx, ...
                                       'COMP4', VEC.C(i,1:4:nx,k)).';
                vx4(i,2:4:nx) = FD_SCH(VEC.V(i,2:4:nx,k), 4*dx, ...
                                       'COMP4', VEC.C(i,2:4:nx,k)).';
                vx4(i,3:4:nx) = FD_SCH(VEC.V(i,3:4:nx,k), 4*dx, ...
                                       'COMP4', VEC.C(i,3:4:nx,k)).';
                vx4(i,4:4:nx) = FD_SCH(VEC.V(i,4:4:nx,k), 4*dx, ...
                                       'COMP4', VEC.C(i,4:4:nx,k)).';
                
                % Compute dw/dx for the full grid (dx).
                wx1(i,:)      = FD_SCH(VEC.W(i,:,k), dx, 'COMP4', ...
                                       VEC.C(i,:,k)).';
                
                % Compute dw/dx for alternative grids skipping 1 node
                % (2*dx).
                wx2(i,1:2:nx) = FD_SCH(VEC.W(i,1:2:nx,k), 2*dx, ...
                                       'COMP4', VEC.C(i,1:2:nx)).';
                wx2(i,2:2:nx) = FD_SCH(VEC.W(i,2:2:nx,k), 2*dx, ...
                                       'COMP4', VEC.C(i,2:2:nx)).';
                
                % Compute dw/dx for alternative grids skipping 3 nodes
                % (4*dx).
                wx4(i,1:4:nx) = FD_SCH(VEC.W(i,1:4:nx,k), 4*dx, ...
                                       'COMP4', VEC.C(i,1:4:nx,k)).';
                wx4(i,2:4:nx) = FD_SCH(VEC.W(i,2:4:nx,k), 4*dx, ...
                                       'COMP4', VEC.C(i,2:4:nx,k)).';
                wx4(i,3:4:nx) = FD_SCH(VEC.W(i,3:4:nx,k), 4*dx, ...
                                       'COMP4', VEC.C(i,3:4:nx,k)).';
                wx4(i,4:4:nx) = FD_SCH(VEC.W(i,4:4:nx,k), 4*dx, ...
                                       'COMP4', VEC.C(i,4:4:nx,k)).';
                
            end
            
            % Compute du/dx, dv/dx and dw/dx using the noise-optimized
            % Richardson scheme.
            DRV.UX(:,:,k) = (1/A)*(A1*ux1 + A2*ux2 + A4*ux4);
            DRV.VX(:,:,k) = (1/A)*(A1*vx1 + A2*vx2 + A4*vx4);
            DRV.WX(:,:,k) = (1/A)*(A1*wx1 + A2*wx2 + A4*wx4);
            
        end
        
        % Initialize the reduced-grid y derivative matrices.
        uy1 = zeros(ny,nx);
        uy2 = zeros(ny,nx);
        uy4 = zeros(ny,nx);
        vy1 = zeros(ny,nx);
        vy2 = zeros(ny,nx);
        vy4 = zeros(ny,nx);
        wy1 = zeros(ny,nx);
        wy2 = zeros(ny,nx);
        wy4 = zeros(ny,nx);
        
        for k = 1:nz
            for j = 1:nx
                
                % Compute du/dy for the full grid (dy).
                uy1(:,j)      = FD_SCH(VEC.U(:,j,k), dy, 'COMP4', ...
                                       VEC.C(:,j,k));
                
                % Compute du/dy for alternative grids skipping 1 node
                % (2*dy).
                uy2(1:2:ny,j) = FD_SCH(VEC.U(1:2:ny,j,k), 2*dy, ...
                                       'COMP4', VEC.C(1:2:ny,j,k));
                uy2(2:2:ny,j) = FD_SCH(VEC.U(2:2:ny,j,k), 2*dy, ...
                                       'COMP4', VEC.C(2:2:ny,j,k));
                
                % Compute du/dy for alternative grids skipping 3 nodes
                % (4*dy).
                uy4(1:4:ny,j) = FD_SCH(VEC.U(1:4:ny,j,k), 4*dy, ...
                                       'COMP4', VEC.C(1:4:ny,j,k));
                uy4(2:4:ny,j) = FD_SCH(VEC.U(2:4:ny,j,k), 4*dy, ...
                                       'COMP4', VEC.C(2:4:ny,j,k));
                uy4(3:4:ny,j) = FD_SCH(VEC.U(3:4:ny,j,k), 4*dy, ...
                                       'COMP4', VEC.C(3:4:ny,j,k));
                uy4(4:4:ny,j) = FD_SCH(VEC.U(4:4:ny,j,k), 4*dy, ...
                                       'COMP4', VEC.C(4:4:ny,j,k));
                
                % Compute dv/dy for the full grid (dy).
                vy1(:,j)      = FD_SCH(VEC.V(:,j,k), dy, 'COMP4', ...
                                       VEC.C(:,j,k));
                
                % Compute dv/dy for alternative grids skipping 1 node
                % (2*dy).
                vy2(1:2:ny,j) = FD_SCH(VEC.V(1:2:ny,j,k), 2*dy, ...
                                       'COMP4', VEC.C(1:2:ny,j,k));
                vy2(2:2:ny,j) = FD_SCH(VEC.V(2:2:ny,j,k), 2*dy, ...
                                       'COMP4', VEC.C(2:2:ny,j,k));
                
                % Compute dv/dy for alternative grids skipping 3 nodes
                % (4*dy).
                vy4(1:4:ny,j) = FD_SCH(VEC.V(1:4:ny,j,k), 4*dy, ...
                                       'COMP4', VEC.C(1:4:ny,j,k));
                vy4(2:4:ny,j) = FD_SCH(VEC.V(2:4:ny,j,k), 4*dy, ...
                                       'COMP4', VEC.C(2:4:ny,j,k));
                vy4(3:4:ny,j) = FD_SCH(VEC.V(3:4:ny,j,k), 4*dy, ...
                                       'COMP4', VEC.C(3:4:ny,j,k));
                vy4(4:4:ny,j) = FD_SCH(VEC.V(4:4:ny,j,k), 4*dy, ...
                                       'COMP4', VEC.C(4:4:ny,j,k));
                
                % Compute dw/dy for the full grid (dy).
                wy1(:,j)      = FD_SCH(VEC.W(:,j,k), dy, 'COMP4', ...
                                       VEC.C(:,j,k));
                
                % Compute dw/dy for alternative grids skipping 1 node
                % (2*dy).
                wy2(1:2:ny,j) = FD_SCH(VEC.W(1:2:ny,j,k), 2*dy, ...
                                       'COMP4', VEC.C(1:2:ny,j,k));
                wy2(2:2:ny,j) = FD_SCH(VEC.W(2:2:ny,j,k), 2*dy, ...
                                       'COMP4', VEC.C(2:2:ny,j,k));
                
                % Compute dw/dy for alternative grids skipping 3 nodes
                % (4*dy).
                wy4(1:4:ny,j) = FD_SCH(VEC.W(1:4:ny,j,k), 4*dy, ...
                                       'COMP4', VEC.C(1:4:ny,j,k));
                wy4(2:4:ny,j) = FD_SCH(VEC.W(2:4:ny,j,k), 4*dy, ...
                                       'COMP4', VEC.C(2:4:ny,j,k));
                wy4(3:4:ny,j) = FD_SCH(VEC.W(3:4:ny,j,k), 4*dy, ...
                                       'COMP4', VEC.C(3:4:ny,j,k));
                wy4(4:4:ny,j) = FD_SCH(VEC.W(4:4:ny,j,k), 4*dy, ...
                                       'COMP4', VEC.C(4:4:ny,j,k));
                
            end
            
            % Compute du/dy, dv/dy and dw/dy using the noise-optimized
            % Richardson scheme.
            DRV.UY(:,:,k) = (1/A)*(A1*uy1 + A2*uy2 + A4*uy4);
            DRV.VY(:,:,k) = (1/A)*(A1*vy1 + A2*vy2 + A4*vy4);
            DRV.WY(:,:,k) = (1/A)*(A1*wy1 + A2*wy2 + A4*wy4);
            
        end
        
        % Initialize the reduced-grid z derivative matrices.
        uz1 = zeros(ny,nx,nz);
        uz2 = zeros(ny,nx,nz);
        uz4 = zeros(ny,nx,nz);
        vz1 = zeros(ny,nx,nz);
        vz2 = zeros(ny,nx,nz);
        vz4 = zeros(ny,nx,nz);
        wz1 = zeros(ny,nx,nz);
        wz2 = zeros(ny,nx,nz);
        wz4 = zeros(ny,nx,nz);
        
        for i = 1:ny
            for j = 1:nx
                
                % Compute du/dz for the full grid (dz).
                uz1(i,j,:)      = FD_SCH(VEC.U(i,j,:), dz, 'COMP4', ...
                                         VEC.C(i,j,:));
                
                % Compute du/dz for alternative grids skipping 1 node
                % (2*dz).
                uz2(i,j,1:2:nz) = FD_SCH(VEC.U(i,j,1:2:nz), 2*dz, ...
                                         'COMP4', VEC.C(i,j,1:2:nz));
                uz2(i,j,2:2:nz) = FD_SCH(VEC.U(i,j,2:2:nz), 2*dz, ...
                                         'COMP4', VEC.C(i,j,2:2:nz));
                
                % Compute du/dz for alternative grids skipping 3 nodes
                % (4*dz).
                uz4(i,j,1:4:nz) = FD_SCH(VEC.U(i,j,1:4:nz), 4*dz, ...
                                         'COMP4', VEC.C(i,j,1:4:nz));
                uz4(i,j,2:4:nz) = FD_SCH(VEC.U(i,j,2:4:nz), 4*dz, ...
                                         'COMP4', VEC.C(i,j,2:4:nz));
                uz4(i,j,3:4:nz) = FD_SCH(VEC.U(i,j,3:4:nz), 4*dz, ...
                                         'COMP4', VEC.C(i,j,3:4:nz));
                uz4(i,j,4:4:nz) = FD_SCH(VEC.U(i,j,4:4:nz), 4*dz, ...
                                         'COMP4', VEC.C(i,j,4:4:nz));
                
                % Compute dv/dz for the full grid (dz).
                vz1(i,j,:)      = FD_SCH(VEC.V(i,j,:), dz, 'COMP4', ...
                                         VEC.C(i,j,:));
                
                % Compute dv/dz for alternative grids skipping 1 node
                % (2*dz).
                vz2(i,j,1:2:nz) = FD_SCH(VEC.V(i,j,1:2:nz), 2*dz, ...
                                         'COMP4', VEC.C(i,j,1:2:nz));
                vz2(i,j,2:2:nz) = FD_SCH(VEC.V(i,j,2:2:nz), 2*dz, ...
                                         'COMP4', VEC.C(i,j,2:2:nz));
                
                % Compute dv/dz for alternative grids skipping 3 nodes
                % (4*dz).
                vz4(i,j,1:4:nz) = FD_SCH(VEC.V(i,j,1:4:nz), 4*dz, ...
                                         'COMP4', VEC.C(i,j,1:4:nz));
                vz4(i,j,2:4:nz) = FD_SCH(VEC.V(i,j,2:4:nz), 4*dz, ...
                                         'COMP4', VEC.C(i,j,2:4:nz));
                vz4(i,j,3:4:nz) = FD_SCH(VEC.V(i,j,3:4:nz), 4*dz, ...
                                         'COMP4', VEC.C(i,j,3:4:nz));
                vz4(i,j,4:4:nz) = FD_SCH(VEC.V(i,j,4:4:nz), 4*dz, ...
                                         'COMP4', VEC.C(i,j,4:4:nz));
                
                % Compute dw/dz for the full grid (dz).
                wz1(i,j,:)      = FD_SCH(VEC.W(i,j,:), dz, 'COMP4', ...
                                         VEC.C(i,j,:));
                
                % Compute dw/dz for alternative grids skipping 1 node
                % (2*dz).
                wz2(i,j,1:2:nz) = FD_SCH(VEC.W(i,j,1:2:nz), 2*dz, ...
                                         'COMP4', VEC.C(i,j,1:2:nz));
                wz2(i,j,2:2:nz) = FD_SCH(VEC.W(i,j,2:2:nz), 2*dz, ...
                                         'COMP4', VEC.C(i,j,2:2:nz));
                
                % Compute dw/dz for alternative grids skipping 3 nodes
                % (4*dz).
                wz4(i,j,1:4:nz) = FD_SCH(VEC.W(i,j,1:4:nz), 4*dz, ...
                                         'COMP4', VEC.C(i,j,1:4:nz));
                wz4(i,j,2:4:nz) = FD_SCH(VEC.W(i,j,2:4:nz), 4*dz, ...
                                         'COMP4', VEC.C(i,j,2:4:nz));
                wz4(i,j,3:4:nz) = FD_SCH(VEC.W(i,j,3:4:nz), 4*dz, ...
                                         'COMP4', VEC.C(i,j,3:4:nz));
                wz4(i,j,4:4:nz) = FD_SCH(VEC.W(i,j,4:4:nz), 4*dz, ...
                                         'COMP4', VEC.C(i,j,4:4:nz));
                
            end
        end
        
        % Compute du/dz, dv/dz and dw/dz using the noise-optimized
        % Richardson scheme.
        DRV.UZ = (1/A)*(A1*uz1 + A2*uz2 + A4*uz4);
        DRV.VZ = (1/A)*(A1*vz1 + A2*vz2 + A4*vz4);
        DRV.WZ = (1/A)*(A1*wz1 + A2*wz2 + A4*wz4);
        
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
