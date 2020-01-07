%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                         ADVECT A GRID OF POINTS                         %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Giuseppe Di Labbio                                                      %
% Department of Mechanical, Industrial & Aerospace Engineering            %
% Concordia University Montréal, Canada                                   %
%                                                                         %
% Last Update: October 5th, 2018 by Giuseppe Di Labbio                    %
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
% A = PP_AdvectGrid(P, VEC, dt);                                          %
% A = PP_AdvectGrid(P, VEC, dt, 'singlestep');                            %
% A = PP_AdvectGrid(P, VEC, dt, 'fracstep', N);                           %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function advects a grid of points with Cartesian coordinates given %
% as a struct (P.X and P.Y) in a velocity field evolving in time. The     %
% user can advect using the raw data (overall time step of 2*dt), a       %
% single time step (dt) or a fraction of a time step. Time-stepping is    %
% performed using the fourth-order Runge-Kutta scheme. This function only %
% applies to two-dimensional data sets for the moment.                    %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'A'          - 1D CELL, ELEMENTS: STRUCTS                               %
%              - One-dimensional cell array of structs each containing    %
%                fields X and Y holding the positions of the advected     %
%                points in time.                                          %
% ----------------------------------------------------------------------- %
% 'dt'         - REAL SCALAR                                              %
%              - Time step of the raw data set.                           %
% ----------------------------------------------------------------------- %
% 'P'          - STRUCT                                                   %
%              - Contains two structs (X and Y) holding the positions of  %
%                the points to be advected.                               %
% ----------------------------------------------------------------------- %
% 'VEC'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - One-dimensional cell array of structs each containing    %
%                information on the spatial mask (C), velocity components %
%                (U and V), and Cartesian grid (X and Y).                 %
% ----------------------------------------------------------------------- %
% Options:                                                                %
% ----------------------------------------------------------------------- %
% 'fracstep'   - SPECIFIC STRING                                          %
%              - Option to use a fractional time step. The number of      %
%                divisions N of the time step must be specified.          %
%              - Default: Not Active                                      %
% ----------------------------------------------------------------------- %
% 'singlestep' - SPECIFIC STRING                                          %
%              - Option to use the full time step.                        %
%              - Default: Not Active                                      %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Advect a grid of points (x,y) = [0.5,1.5]x[0.25:0.75] with a constant   %
% grid spacing of 0.1 in a time-dependent double gyre over the domain     %
% (x,y) = [0,2]x[0,1] with a constant grid spacing of 0.01 over the time  %
% interval [0,20] with time-step size 0.1. Use A = 0.1, epsilon = 0.25,   %
% and omega = 2*pi/10. Use the time-step dt = 0.1 for the advection (i.e. %
% use the 'singlestep' option).                                           %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A    = 0.1;                                                          %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT]    = GEN_DoubleGyre(x, y, t, A, epsn, omga);              %
% >> [grd.X,grd.Y] = meshgrid(0.5:0.1:1.5,0.25:0.1:0.75);                 %
% >> A = PP_AdvectGrid(grd, VEC, 0.1, 'singlestep');                      %
% >> for k = 1:length(t)                                                  %
% quiver(VEC{k}.X(1:4:end,1:4:end), VEC{k}.Y(1:4:end,1:4:end), ...        %
%        VEC{k}.U(1:4:end,1:4:end), VEC{k}.V(1:4:end,1:4:end));           %
% axis([0 2 0 1]);                                                        %
% hold on;                                                                %
% scatter(mat2vec(A{k}.X), mat2vec(A{k}.Y), 'r', 'filled');               %
% hold off;                                                               %
% pause(0.25);                                                            %
% end                                                                     %
% >> clear k;                                                             %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Advect a circular grid of points with center (1.5,0.5), a radius of 0.1 %
% and a grid spacing of 0.001 in a time-dependent double gyre over the    %
% domain (x,y) = [0,2]x[0,1] with a constant grid spacing of 0.01 over    %
% the time interval [0,20] with time-step size 0.1. Use A = 0.1, epsilon  %
% = 0.25, and omega = 2*pi/10. Use the time-step dt = 0.1 for the         %
% advection (i.e. use the 'singlestep' option).                           %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A    = 0.1;                                                          %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> grd = PP_CreateBlobs([1.0 0.5], 0.1, 0.001);                         %
% >> A = PP_AdvectGrid(grd, VEC, 0.1, 'singlestep');                      %
% >> for k = 1:length(t)                                                  %
% quiver(VEC{k}.X(1:4:end,1:4:end), VEC{k}.Y(1:4:end,1:4:end), ...        %
%        VEC{k}.U(1:4:end,1:4:end), VEC{k}.V(1:4:end,1:4:end));           %
% axis([0 2 0 1]);                                                        %
% hold on;                                                                %
% scatter(mat2vec(A{k}.X), mat2vec(A{k}.Y), 'r', 'filled');               %
% hold off;                                                               %
% pause(0.25);                                                            %
% end                                                                     %
% >> clear k;                                                             %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% N/A                                                                     %
%                                                                         %
% Called in:                                                              %
% PP_FTLE                                                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_AdvectGrid
function [A] = PP_AdvectGrid(P, VEC, dt, varargin)

% Determine the number of time steps over which to advect particles.
n = length(VEC);

% Determine the number of particles to advect in the x and y directions.
[ny,nx] = size(P.X);

if nargin == 4 && strcmpi(varargin{1}, 'singlestep')
    
    % Initialize a cell array of structs (A) containing the X and Y
    % coordinates of the advected grid in time.
    A = cell(n,1);
    A{1} = struct('X', P.X, 'Y', P.Y);
    for k = 2:length(A)
        A{k} = struct('X', zeros(ny,nx), 'Y', zeros(ny,nx));
    end
    
    c = 1;
    for k = 1:n-1
        
        Xhalf = 0.5*(VEC{k}.X + VEC{k+1}.X);
        Yhalf = 0.5*(VEC{k}.Y + VEC{k+1}.Y);
        Uhalf = 0.5*(VEC{k}.U + VEC{k+1}.U);
        Vhalf = 0.5*(VEC{k}.V + VEC{k+1}.V);
        
        K1u = interp2(VEC{k}.X, VEC{k}.Y, VEC{k}.U, A{c}.X, A{c}.Y, ...
                      'cubic');
        K1v = interp2(VEC{k}.X, VEC{k}.Y, VEC{k}.V, A{c}.X, A{c}.Y, ...
                      'cubic');
        K2u = interp2(Xhalf, Yhalf, Uhalf, A{c}.X + 0.5*dt*K1u, ...
                      A{c}.Y + 0.5*dt*K1v, 'cubic');
        K2v = interp2(Xhalf, Yhalf, Vhalf, A{c}.X + 0.5*dt*K1u, ...
                      A{c}.Y + 0.5*dt*K1v, 'cubic');
        K3u = interp2(Xhalf, Yhalf, Uhalf, A{c}.X + 0.5*dt*K2u, ...
                      A{c}.Y + 0.5*dt*K2v, 'cubic');
        K3v = interp2(Xhalf, Yhalf, Vhalf, A{c}.X + 0.5*dt*K2u, ...
                      A{c}.Y + 0.5*dt*K2v, 'cubic');
        K4u = interp2(VEC{k+1}.X, VEC{k+1}.Y, VEC{k+1}.U, ...
                      A{c}.X + dt*K3u, A{c}.Y + dt*K3v, 'cubic');
        K4v = interp2(VEC{k+1}.X, VEC{k+1}.Y, VEC{k+1}.V, ...
                      A{c}.X + dt*K3u, A{c}.Y + dt*K3v, 'cubic');
        
        A{c+1}.X  = A{c}.X + (dt/6)*(K1u + 2*K2u + 2*K3u + K4u);
        A{c+1}.Y  = A{c}.Y + (dt/6)*(K1v + 2*K2v + 2*K3v + K4v);
        c = c + 1;
        
    end
    
elseif nargin == 3
   
    % Initialize a cell array of structs (A) containing the X and Y
    % coordinates of the advected grid in time.
    A = cell(ceil(n/2),1);
    A{1} = struct('X', P.X, 'Y', P.Y);
    for k = 2:length(A)
        A{k} = struct('X', zeros(ny,nx), 'Y', zeros(ny,nx));
    end
    
    c = 1;
    for k = 1:2:n-2
        
        K1u = interp2(VEC{k}.X, VEC{k}.Y, VEC{k}.U, A{c}.X, A{c}.Y, ...
                      'cubic');
        K1v = interp2(VEC{k}.X, VEC{k}.Y, VEC{k}.V, A{c}.X, A{c}.Y, ...
                      'cubic');
        K2u = interp2(VEC{k+1}.X, VEC{k+1}.Y, VEC{k+1}.U, ...
                      A{c}.X + dt*K1u, A{c}.Y + dt*K1v, 'cubic');
        K2v = interp2(VEC{k+1}.X, VEC{k+1}.Y, VEC{k+1}.V, ...
                      A{c}.X + dt*K1u, A{c}.Y + dt*K1v, 'cubic');
        K3u = interp2(VEC{k+1}.X, VEC{k+1}.Y, VEC{k+1}.U, ...
                      A{c}.X + dt*K2u, A{c}.Y + dt*K2v, 'cubic');
        K3v = interp2(VEC{k+1}.X, VEC{k+1}.Y, VEC{k+1}.V, ...
                      A{c}.X + dt*K2u, A{c}.Y + dt*K2v, 'cubic');
        K4u = interp2(VEC{k+2}.X, VEC{k+2}.Y, VEC{k+2}.U, ...
                      A{c}.X + 2*dt*K3u, A{c}.Y + 2*dt*K3v, 'cubic');
        K4v = interp2(VEC{k+2}.X, VEC{k+2}.Y, VEC{k+2}.V, ...
                      A{c}.X + 2*dt*K3u, A{c}.Y + 2*dt*K3v, 'cubic');
        
        A{c+1}.X  = A{c}.X + (2*dt/6)*(K1u + 2*K2u + 2*K3u + K4u);
        A{c+1}.Y  = A{c}.Y + (2*dt/6)*(K1v + 2*K2v + 2*K3v + K4v);
        c = c + 1;
        
    end
    
elseif nargin == 5 && strcmpi(varargin{1}, 'fracstep')
    
    N = varargin{2};
    
    % Initialize a cell array of structs (A) containing the X and Y
    % coordinates of the advected grid in time.
    A = cell(N*n,1);
    A{1} = struct('X', P.X, 'Y', P.Y);
    for k = 2:length(A)
        A{k} = struct('X', zeros(ny,nx), 'Y', zeros(ny,nx));
    end
    
    c = 1;
    for k1 = 1:n-1
        
        X = permute(cat(3,VEC{k1}.X,VEC{k1+1}.X),[3 1 2]);
        Y = permute(cat(3,VEC{k1}.Y,VEC{k1+1}.Y),[3 1 2]);
        U = permute(cat(3,VEC{k1}.U,VEC{k1+1}.U),[3 1 2]);
        V = permute(cat(3,VEC{k1}.V,VEC{k1+1}.V),[3 1 2]);
        Xfrac = permute(interp1([0 1], X, 0:1/(2*N):1), [2 3 1]);
        Yfrac = permute(interp1([0 1], Y, 0:1/(2*N):1), [2 3 1]);
        Ufrac = permute(interp1([0 1], U, 0:1/(2*N):1), [2 3 1]);
        Vfrac = permute(interp1([0 1], V, 0:1/(2*N):1), [2 3 1]);
        for k2 = 1:2:2*N
            
            K1u = interp2(Xfrac(:,:,k2), Yfrac(:,:,k2), Ufrac(:,:,k2),  ...
                          A{c}.X, A{c}.Y, 'cubic');
            K1v = interp2(Xfrac(:,:,k2), Yfrac(:,:,k2), Vfrac(:,:,k2),  ...
                          A{c}.X, A{c}.Y, 'cubic');
            K2u = interp2(Xfrac(:,:,k2+1), Yfrac(:,:,k2+1),             ...
                          Ufrac(:,:,k2+1), A{c}.X + 0.5*(dt/N)*K1u,     ...
                          A{c}.Y + 0.5*(dt/N)*K1v, 'cubic');
            K2v = interp2(Xfrac(:,:,k2+1), Yfrac(:,:,k2+1),             ...
                          Vfrac(:,:,k2+1), A{c}.X + 0.5*(dt/N)*K1u,     ...
                          A{c}.Y + 0.5*(dt/N)*K1v, 'cubic');
            K3u = interp2(Xfrac(:,:,k2+1), Yfrac(:,:,k2+1),             ...
                          Ufrac(:,:,k2+1), A{c}.X + 0.5*(dt/N)*K2u,     ...
                          A{c}.Y + 0.5*(dt/N)*K2v, 'cubic');
            K3v = interp2(Xfrac(:,:,k2+1), Yfrac(:,:,k2+1),             ...
                          Vfrac(:,:,k2+1), A{c}.X + 0.5*(dt/N)*K2u,     ...
                          A{c}.Y + 0.5*(dt/N)*K2v, 'cubic');
            K4u = interp2(Xfrac(:,:,k2+2), Yfrac(:,:,k2+2), ...
                          Ufrac(:,:,k2+2), A{c}.X + (dt/N)*K3u, ...
                          A{c}.Y + (dt/N)*K3v, 'cubic');
            K4v = interp2(Xfrac(:,:,k2+2), Yfrac(:,:,k2+2), ...
                          Vfrac(:,:,k2+2), A{c}.X + (dt/N)*K3u, ...
                          A{c}.Y + (dt/N)*K3v, 'cubic');
                  
            A{c+1}.X  = A{c}.X + (dt/6)*(K1u + 2*K2u + 2*K3u + K4u)/N;
            A{c+1}.Y  = A{c}.Y + (dt/6)*(K1v + 2*K2v + 2*K3v + K4v)/N;
            c = c + 1;
        
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
