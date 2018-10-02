%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                      GAMMA-2 VORTEX IDENTIFICATION                      %
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
% GAM2 = PP_Gamma2(VEC, S, rf);                                           %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function identifies a vortex center using the Gamma-2 method       %
% defined by Michard & Graftieaux (2001) (see Graftieaux, L., Michard, M. %
% & Grosjean, N. (2001). Combining PIV, POD and vortex identification     %
% algorithms for the study of unsteady turbulent swirling flows.          %
% Measurement Science and Technology, 12, 1422-1429). Note that Gamma-2   %
% is Galilean invariant whereas Gamma-1 is not.                           %
%                                                                         %
% References:                                                             %
% Coming soon ...                                                         %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'GAM2'       - STRUCT                                                   %
%              - Gamma-2 struct containg the Cartesian coordinates 'X'    %
%                'Y' as well as the value of Gamma-1 in 'Val'.            %
% ----------------------------------------------------------------------- %
% 'rf'         - INTEGER SCALAR                                           %
%              - The refinement factor denotes the multiple over which    %
%                grid is to be refined.                                   %
% ----------------------------------------------------------------------- %
% 'S'          - INTEGER SCALAR                                           %
%              - Window size to use in calculating Gamma-2.               %
% ----------------------------------------------------------------------- %
% 'VEC'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - One-dimensional cell array of structs each containing    %
%                information on the spatial mask (C), velocity components %
%                (U and V), and Cartesian grid (X and Y).                 %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Track the two vortex centers for a time-dependent double gyre on the    %
% domain (x,y) = [0,2]x[0,1] with a constant grid spacing of 0.01 over    %
% the time interval [0,20] with time-step size 0.1. Use A = 0.1, epsilon  %
% = 0.25, and omega = 2*pi/10. Use a refinement factor and window size    %
% equal to 4.                                                             %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> GAM2 = cell(length(t),1);                                            %
% >> P1 = zeros(length(t),2);                                             %
% >> P2 = zeros(length(t),2);                                             %
% >> for k = 1:length(t)                                                  %
% GAM2{k} = PP_Gamma2(VEC{k}, 4, 4);                                      %
% [Mtemp,Mind1] = min(GAM2{k}.Val);                                       %
% [~, Mind2] = min(Mtemp);                                                %
% P1(k,1) = GAM2{k}.X(Mind1(Mind2), Mind2);                               %
% P1(k,2) = GAM2{k}.Y(Mind1(Mind2), Mind2);                               %
% [Mtemp, Mind1] = max(GAM2{k}.Val);                                      %
% [~, Mind2] = max(Mtemp);                                                %
% P2(k,1) = GAM2{k}.X(Mind1(Mind2), Mind2);                               %
% P2(k,2) = GAM2{k}.Y(Mind1(Mind2), Mind2);                               %
% end                                                                     %
% >> clear Mind1 Mind2 Mtemp;                                             %
% >> for k = 1:length(t)                                                  %
% quiver(VEC{k}.X(1:4:end,1:4:end), VEC{k}.Y(1:4:end,1:4:end), ...        %
%        VEC{k}.U(1:4:end,1:4:end), VEC{k}.V(1:4:end,1:4:end), 'k');      %
% hold on;                                                                %
% scatter(P1(k,1), P1(k,2), 'b', 'filled');                               %
% scatter(P2(k,1), P2(k,2), 'r', 'filled');                               %
% hold off;                                                               %
% axis([0 2 0 1]);                                                        %
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
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_Gamma2
function [GAM2] = PP_Gamma2(VEC, S, rf)

% Determine the new grid size in x and y.
nx = rf*size(VEC.X, 2);
ny = rf*size(VEC.Y, 1);

% Define a refined grid to compute the Gamma-1 field over.
x = linspace(VEC.X(1,1), VEC.X(1,end), nx);
y = fliplr(linspace(VEC.Y(end,1), VEC.Y(1,1), ny));
[P.X, P.Y] = meshgrid(x, y);

U = interp2(VEC.X, VEC.Y, VEC.U, P.X, P.Y, 'cubic');
V = interp2(VEC.X, VEC.Y, VEC.V, P.X, P.Y, 'cubic');

Nx = floor(nx/S);
Ny = floor(ny/S);
if nx/S == Nx && ny/S == Ny
    GAM2 = struct('X', zeros(Ny, Nx), 'Y', zeros(Ny, Nx), ...
                  'Val', zeros(Ny, Nx));
    flag = 0;
elseif nx/S == Nx && ny/S ~= Ny
    GAM2 = struct('X', zeros(Ny+1, Nx), 'Y', zeros(Ny+1, Nx), ...
                  'Val', zeros(Ny+1, Nx));
    flag = 1;
elseif nx/S ~= Nx && ny/S == Ny
    GAM2 = struct('X', zeros(Ny, Nx+1), 'Y', zeros(Ny, Nx+1), ...
                  'Val', zeros(Ny, Nx+1));
    flag = 2;
else
    GAM2 = struct('X', zeros(Ny+1, Nx+1), 'Y', zeros(Ny+1, Nx+1), ...
                  'Val', zeros(Ny+1, Nx+1));
    flag = 3;
end

r = [0 0 0];
um = [0 0 0];
up = [0 0 0];
for i = 1:Ny
    for j = 1:Nx
        
        up(1) = trapz(trapz(U(S*(i-1)+1:S*i,S*(j-1)+1:S*j),1),2)/S^2;
        up(2) = trapz(trapz(V(S*(i-1)+1:S*i,S*(j-1)+1:S*j),1),2)/S^2;
        GAM2.X(i,j) = (x(S*(j-1)+1) + x(S*j))/2;
        GAM2.Y(i,j) = (y(S*(i-1)+1) + y(S*i))/2;
        for ii = 1:S
            for jj = 1:S
                
                r(1)  = P.X(S*(i-1)+ii,S*(j-1)+jj) - GAM2.X(i,j);
                r(2)  = P.Y(S*(i-1)+ii,S*(j-1)+jj) - GAM2.Y(i,j);
                um(1) = U(S*(i-1)+ii,S*(j-1)+jj) - up(1);
                um(2) = V(S*(i-1)+ii,S*(j-1)+jj) - up(2);
                C     = dot(cross(r,um), [0 0 1])/(norm(um)*norm(r));
                GAM2.Val(i,j) = GAM2.Val(i,j) + C;
                
            end
        end
        GAM2.Val(i,j) = GAM2.Val(i,j)/S^2;
        
    end
end

switch flag
    case 1
        for j = 1:Nx
            up(1) = trapz(trapz(U(S*Ny+1:end,S*(j-1)+1:S*j),1),2)/...
                         (S*(ny-S*Ny));
            up(2) = trapz(trapz(V(S*Ny+1:end,S*(j-1)+1:S*j),1),2)/...
                         (S*(ny-S*Ny));
            GAM2.X(Ny+1,j) = (x(S*(j-1)+1) + x(S*j))/2;
            GAM2.Y(Ny+1,j) = (y(S*Ny+1) + y(end))/2;
            for ii = 1:ny-S*Ny
                for jj = 1:S
                    
                    r(1)  = P.X(S*Ny+ii,S*(j-1)+jj) - GAM2.X(Ny+1,j);
                    r(2)  = P.Y(S*Ny+ii,S*(j-1)+jj) - GAM2.Y(Ny+1,j);
                    um(1) = U(S*Ny+ii,S*(j-1)+jj) - up(1);
                    um(2) = V(S*Ny+ii,S*(j-1)+jj) - up(2);
                    C     = dot(cross(r,um), [0 0 1])/(norm(um)*norm(r));
                    GAM2.Val(Ny+1,j) = GAM2.Val(Ny+1,j) + C;
                    
                end
            end
            GAM2.Val(Ny+1,j) = GAM2.Val(Ny+1,j)/(S*(ny-S*Ny));
        end
        
    case 2
        for i = 1:Ny
            up(1) = trapz(trapz(U(S*(i-1)+1:S*i,S*Nx+1:end),1),2)/...
                         (S*(nx-S*Nx));
            up(2) = trapz(trapz(V(S*(i-1)+1:S*i,S*Nx+1:end),1),2)/...
                         (S*(nx-S*Nx));
            GAM2.X(i,Nx+1) = (x(S*Nx+1) + x(end))/2;
            GAM2.Y(i,Nx+1) = (y(S*(i-1)+1) + y(S*i))/2;
            for ii = 1:S
                for jj = 1:nx-S*Nx
                    
                    r(1)  = P.X(S*(i-1)+ii,S*Nx+jj) - GAM2.X(i,Nx+1);
                    r(2)  = P.Y(S*(i-1)+ii,S*Nx+jj) - GAM2.Y(i,Nx+1);
                    um(1) = U(S*(i-1)+ii,S*Nx+jj) - up(1);
                    um(2) = V(S*(i-1)+ii,S*Nx+jj) - up(2);
                    C     = dot(cross(r,um), [0 0 1])/(norm(um)*norm(r));
                    GAM2.Val(i,Nx+1) = GAM2.Val(i,Nx+1) + C;
                    
                end
            end
            GAM2.Val(i,Nx+1) = GAM2.Val(i,Nx+1)/(S*(nx-S*Nx));
        end
        
    case 3
        for j = 1:Nx
            up(1) = trapz(trapz(U(S*Ny+1:end,S*(j-1)+1:S*j),1),2)/...
                         (S*(ny-S*Ny));
            up(2) = trapz(trapz(V(S*Ny+1:end,S*(j-1)+1:S*j),1),2)/...
                         (S*(ny-S*Ny));
            GAM2.X(Ny+1,j) = (x(S*(j-1)+1) + x(S*j))/2;
            GAM2.Y(Ny+1,j) = (y(S*Ny+1) + y(end))/2;
            for ii = 1:ny-S*Ny
                for jj = 1:S
                    
                    r(1)  = P.X(S*Ny+ii,S*(j-1)+jj) - GAM2.X(Ny+1,j);
                    r(2)  = P.Y(S*Ny+ii,S*(j-1)+jj) - GAM2.Y(Ny+1,j);
                    um(1) = U(S*Ny+ii,S*(j-1)+jj) - up(1);
                    um(2) = V(S*Ny+ii,S*(j-1)+jj) - up(2);
                    C     = dot(cross(r,um), [0 0 1])/(norm(um)*norm(r));
                    GAM2.Val(Ny+1,j) = GAM2.Val(Ny+1,j) + C;
                    
                end
            end
            GAM2.Val(Ny+1,j) = GAM2.Val(Ny+1,j)/(S*(ny-S*Ny));
        end
        
        for i = 1:Ny
            up(1) = trapz(trapz(U(S*(i-1)+1:S*i,S*Nx+1:end),1),2)/...
                         (S*(nx-S*Nx));
            up(2) = trapz(trapz(V(S*(i-1)+1:S*i,S*Nx+1:end),1),2)/...
                         (S*(nx-S*Nx));
            GAM2.X(i,Nx+1) = (x(S*Nx+1) + x(end))/2;
            GAM2.Y(i,Nx+1) = (y(S*(i-1)+1) + y(S*i))/2;
            for ii = 1:S
                for jj = 1:nx-S*Nx
                    
                    r(1)  = P.X(S*(i-1)+ii,S*Nx+jj) - GAM2.X(i,Nx+1);
                    r(2)  = P.Y(S*(i-1)+ii,S*Nx+jj) - GAM2.Y(i,Nx+1);
                    um(1) = U(S*(i-1)+ii,S*Nx+jj) - up(1);
                    um(2) = V(S*(i-1)+ii,S*Nx+jj) - up(2);
                    C     = dot(cross(r,um), [0 0 1])/(norm(um)*norm(r));
                    GAM2.Val(i,Nx+1) = GAM2.Val(i,Nx+1) + C;
                    
                end
            end
            GAM2.Val(i,Nx+1) = GAM2.Val(i,Nx+1)/(S*(nx-S*Nx));
        end
        
        up(1) = trapz(trapz(U(S*Ny+1:end,S*(j-1)+1:S*j),1),2)/...
                     ((ny-S*Ny)*(nx-S*Nx));
        up(2) = trapz(trapz(V(S*Ny+1:end,S*(j-1)+1:S*j),1),2)/...
                     ((ny-S*Ny)*(nx-S*Nx));
        GAM2.X(Ny+1,Nx+1) = (x(S*Nx+1) + x(end))/2;
        GAM2.Y(Ny+1,Nx+1) = (y(S*Ny+1) + y(end))/2;
        for ii = 1:ny-S*Ny
            for jj = 1:nx-S*Nx
                
                r(1)  = P.X(S*Ny+ii,S*Nx+jj) - GAM2.X(Ny+1,Nx+1);
                r(2)  = P.Y(S*Ny+ii,S*Nx+jj) - GAM2.Y(Ny+1,Nx+1);
                um(1) = U(S*Ny+ii,S*Nx+jj) - up(1);
                um(2) = V(S*Ny+ii,S*Nx+jj) - up(2);
                C     = dot(cross(r,um), [0 0 1])/(norm(um)*norm(r));
                GAM2.Val(Ny+1,Nx+1) = GAM2.Val(Ny+1,Nx+1) + C;
                
            end
        end
        GAM2.Val(Ny+1,Nx+1) = GAM2.Val(Ny+1,Nx+1)/((ny-S*Ny)*(nx-S*Nx));
        
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
