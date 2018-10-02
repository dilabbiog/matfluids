%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                 ADVECT SELECT POINTS (COMPUTE PATHLINES)                %
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
% A = PP_AdvectPoints(P, VEC, dt);                                        %
% A = PP_AdvectPoints(P, VEC, dt, 'singlestep');                          %
% A = PP_AdvectPoints(P, VEC, dt, 'fracstep', frac);  - N/A -             %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function advects a list of points with Cartesian coordinates given %
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
% 'A'          - STRUCT                                                   %
%              - Contains two structs (X and Y) holding the positions of  %
%                the initial points in the first row and their time       %
%                series in their respective columns.                      %
% ----------------------------------------------------------------------- %
% 'dt'         - REAL SCALAR                                              %
%              - Time step of the raw data set.                           %
% ----------------------------------------------------------------------- %
% 'P'          - STRUCT                                                   %
%              - Contains two structs (X and Y) holding the positions of  %
%                the points to be advected as a list of points (column or %
%                row vectors will both work.                              %
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
% EXAMPLE                                                                 %
%                                                                         %
% Advect the points (x,y) = (0,-0.015), (0,-0.010), (0,-0.005), (0,0),    %
% (0,0.005), (0,0.010) and (0,0.015) in a steady Hagen-Poiseuille flow    %
% in a circular pipe of length 10 m and radius 2 cm. Use a grid spacing   %
% of 0.02 cm in the radial direction and 0.1 m in the axial direction.    %
% Assume a dynamic viscosity of 1 cP and a pressure difference of 0.1 kPa %
% across the length of the pipe. Use a time step of dt = 0.1 s to advect  %
% the particles from t = 0 s to t = 5 s (use the 'singlestep' option).    %
%                                                                         %
% >> R  = 0.02;                                                           %
% >> dr = 0.0002;                                                         %
% >> L  = 10;                                                             %
% >> dl = 0.1;                                                            %
% >> t  = (0:0.1:5).';                                                    %
% >> mu = 0.001;                                                          %
% >> dP = 100;                                                            %
% >> [VEC, VGT] = GEN_HagenPoiseuille(R, dr, L, dl, t, mu, dP);           %
% >> P.X = [0 0 0 0 0 0 0];                                               %
% >> P.Y = [-0.015 -0.010 -0.005 0 0.005 0.010 0.015];                    %
% >> A = PP_AdvectPoints(P, VEC, 0.1, 'singlestep');                      %
% >> for k = 1:length(t)                                                  %
% quiver(VEC{1}.X(1:4:end,1:4:end), VEC{1}.Y(1:4:end,1:4:end), ...        %
%        VEC{1}.U(1:4:end,1:4:end), VEC{1}.V(1:4:end,1:4:end), ...        %
%        'ShowArrowHead', 'off');                                         %
% axis([0 10 -0.05 0.05]);                                                %
% hold on;                                                                %
% plot([0 10], [0.02 0.02],   'LineWidth', 3, 'Color', 'k');              %
% plot([0 10], [-0.02 -0.02], 'LineWidth', 3, 'Color', 'k');              %
% scatter(A.X(k,:), A.Y(k,:), 'k');                                       %
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
% PP_Streaklines                                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_AdvectPoints
function [A] = PP_AdvectPoints(P, VEC, dt, varargin)

% Determine the number of time steps over which to advect particles.
n = length(VEC);

% Determine the number of particles to advect.
p = length(P.X);

if nargin == 4 && strcmpi(varargin{1},'singlestep')
    
    % Initialize a struct (A) containing the X and Y coordinates of the
    % advected points in time.
    A = struct('X', zeros(n, p), 'Y', zeros(n, p));
    
    c = 1;
    A.X(1,:) = P.X;
    A.Y(1,:) = P.Y;
    for k = 1:n-1
        
        Xhalf = 0.5*(VEC{k}.X + VEC{k+1}.X);
        Yhalf = 0.5*(VEC{k}.Y + VEC{k+1}.Y);
        Uhalf = 0.5*(VEC{k}.U + VEC{k+1}.U);
        Vhalf = 0.5*(VEC{k}.V + VEC{k+1}.V);
        
        K1u = interp2(VEC{k}.X, VEC{k}.Y, VEC{k}.U, A.X(c,:), A.Y(c,:), ...
                      'cubic');
        K1v = interp2(VEC{k}.X, VEC{k}.Y, VEC{k}.V, A.X(c,:), A.Y(c,:), ...
                      'cubic');
        K2u = interp2(Xhalf, Yhalf, Uhalf, A.X(c,:) + 0.5*dt*K1u, ...
                      A.Y(c,:) + 0.5*dt*K1v, 'cubic');
        K2v = interp2(Xhalf, Yhalf, Vhalf, A.X(c,:) + 0.5*dt*K1u, ...
                      A.Y(c,:) + 0.5*dt*K1v, 'cubic');
        K3u = interp2(Xhalf, Yhalf, Uhalf, A.X(c,:) + 0.5*dt*K2u, ...
                      A.Y(c,:) + 0.5*dt*K2v, 'cubic');
        K3v = interp2(Xhalf, Yhalf, Vhalf, A.X(c,:) + 0.5*dt*K2u, ...
                      A.Y(c,:) + 0.5*dt*K2v, 'cubic');
        K4u = interp2(VEC{k+1}.X, VEC{k+1}.Y, VEC{k+1}.U, ...
                      A.X(c,:) + dt*K3u, A.Y(c,:) + dt*K3v, 'cubic');
        K4v = interp2(VEC{k+1}.X, VEC{k+1}.Y, VEC{k+1}.V, ...
                      A.X(c,:) + dt*K3u, A.Y(c,:) + dt*K3v, 'cubic');
        
        A.X(c+1,:)  = A.X(c,:) + (dt/6)*(K1u + 2*K2u + 2*K3u + K4u);
        A.Y(c+1,:)  = A.Y(c,:) + (dt/6)*(K1v + 2*K2v + 2*K3v + K4v);
        c = c + 1;
        
    end
    
else

    % Initialize a struct (A) containing the X and Y coordinates of the
    % advected points in time.
    A = struct('X', zeros(ceil(n/2), p), 'Y', zeros(ceil(n/2), p));
    
    c = 1;
    A.X(1,:) = P.X;
    A.Y(1,:) = P.Y;
    for k = 1:2:n-2
        
        K1u = interp2(VEC{k}.X, VEC{k}.Y, VEC{k}.U, A.X(c,:), A.Y(c,:), ...
                      'cubic');
        K1v = interp2(VEC{k}.X, VEC{k}.Y, VEC{k}.V, A.X(c,:), A.Y(c,:), ...
                      'cubic');
        K2u = interp2(VEC{k+1}.X, VEC{k+1}.Y, VEC{k+1}.U, ...
                      A.X(c,:) + dt*K1u, A.Y(c,:) + dt*K1v, 'cubic');
        K2v = interp2(VEC{k+1}.X, VEC{k+1}.Y, VEC{k+1}.V, ...
                      A.X(c,:) + dt*K1u, A.Y(c,:) + dt*K1v, 'cubic');
        K3u = interp2(VEC{k+1}.X, VEC{k+1}.Y, VEC{k+1}.U, ...
                      A.X(c,:) + dt*K2u, A.Y(c,:) + dt*K2v, 'cubic');
        K3v = interp2(VEC{k+1}.X, VEC{k+1}.Y, VEC{k+1}.V, ...
                      A.X(c,:) + dt*K2u, A.Y(c,:) + dt*K2v, 'cubic');
        K4u = interp2(VEC{k+2}.X, VEC{k+2}.Y, VEC{k+2}.U, ...
                      A.X(c,:) + 2*dt*K3u, A.Y(c,:) + 2*dt*K3v, 'cubic');
        K4v = interp2(VEC{k+2}.X, VEC{k+2}.Y, VEC{k+2}.V, ...
                      A.X(c,:) + 2*dt*K3u, A.Y(c,:) + 2*dt*K3v, 'cubic');
        
        A.X(c+1,:)  = A.X(c,:) + (2*dt/6)*(K1u + 2*K2u + 2*K3u + K4u);
        A.Y(c+1,:)  = A.Y(c,:) + (2*dt/6)*(K1v + 2*K2v + 2*K3v + K4v);
        c = c + 1;
        
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
