%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       MODAL DECOMPOSITION TOOLBOX                       %
%                                                                         %
% podSpaceOnly                                                            %
% Proper Orthogonal Decomposition                                         %
% Space-only proper orthogonal decomposition method                       %
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
% Copyright (C) 2023 Giuseppe Di Labbio                                   %
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
% [Phi, A, props] = podSpaceOnly(X);                                      %
% [Phi, A, props] = podSpaceOnly(X, 'zeromean');                          %
% [Phi, A, props] = podSpaceOnly(___, Name, Value);                       %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Compute the space-only proper orthogonal decomposition [1-3] of a       %
% dataset. The decomposition uses the common "method of snapshots" [1].   %
% The input data is stored in a matrix X whose columns represent distinct %
% realizations of a measurement (e.g., snapshots over time) and whose     %
% rows represent a series of measurements (e.g., taken at a given time).  %
% The decomposition can be computed at a reduced rank by specifying a     %
% tolerance (a cutoff threshold based on the fraction of the largest      %
% eigenvalue).                                                            %
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
% [1] Sirovich, L. (1987). Turbulence and the dynamics of coherent        %
%     structures. I. Coherent structures. Quarterly of Applied            %
%     Mathematics, 45(3), 561-571.                                        %
% [2] Lumley, J. L. (1967). The structure of inhomogeneous turbulent      %
%     flows. In Atmospheric Turbulence and Radio Propagation (ed. A. M.   %
%     Yaglom & V. I. Tatarski), pp. 166-178. Nauka.                       %
% [3] Towne, A., Schmidt, O. T., & Colonius, T. (2018). Spectral proper   %
%     orthogonal decomposition and its relationship to dynamic mode       %
%     decomposition and resolvent analysis. Journal of Fluid Mechanics,   %
%     847, 821-867.
%                                                                         %
% ======================================================================= %
% Input Arguments (Required):                                             %
% ----------------------------------------------------------------------- %
% 'X'            LOGICAL/NUMERIC 2-DIMENSIONAL ARRAY                      %
%              ~ Dataset on which to compute the space-only proper        %
%                orthogonal decomposition. In the case of fluid flow, the %
%                columns of X often represent sequential snapshots while  %
%                the rows often represent the velocity components at      %
%                various points in the flow (within a given snapshot).    %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'zeromean'     POSITIVE REAL SCALAR                                     %
%                Default: 0                                               %
%              ~ Subtract the mean of all realizations from the dataset.  %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% 'rankTol'      NONNEGATIVE REAL SCALAR                                  %
%                Default: 0                                               %
%              ~ Tolerance to reduce the rank of the algorithm at the     %
%                eigenvalue problem level. The tolerance is taken as a    %
%                fraction of the largest eigenvalue and is used as a      %
%                cutoff (EV < rankTol*max(EV) are ignored).               %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'Phi'          COMPLEX 2-DIMENSIONAL ARRAY                              %
%              ~ Space-only proper orthogonal modes of the dataset. The   %
%                number of rows is identical to that of the input data    %
%                whereas the number of columns depends on rank reduction. %
% ----------------------------------------------------------------------- %
% 'A'            COMPLEX 2-DIMENSIONAL ARRAY                              %
%              ~ Time dynamics of the space-only proper orthogonal modes  %
%                of the dataset. Each row represents the time dynamics of %
%                the mode in the respective column of Phi. The number of  %
%                columns is therefore identical to that of the input data %
%                (i.e., the number of realizations).                      %
% ----------------------------------------------------------------------- %
% 'props'        STRUCT ARRAY (1 X 1)                                     %
%              ~ Properties of the space-only proper orthogonal mode      %
%                decomposition. The structure array contains the          %
%                following fields:                                        %
%                 1) E  : Mode energy                                     %
%                 2) H  : Shannon entropy of modal energy distribution    %
%                 3) rk : Rank of the space-only POD                      %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Consider a double gyre on the domain (x,y) = ([0,2],[0,1]) with a       %
% constant grid spacing of 0.01 over the time interval [0,20] with time   %
% step size 0.1 (A = 0.1, epsi = 0.25, omega = 2*pi/10). Compute the      %
% space-only proper orthogonal decomposition of the velocity field. Plot  %
% the first spatial mode (x component) and its time dynamics. Plot the    %
% energy spectrum for the first 20 modes, normalized by the energy of the %
% first mode.                                                             %
%                                                                         %
% >> N       = struct('x', 201, 'y', 101, 't', 201);                      %
% >> coord.x = linspace(0,  2, N.x).';                                    %
% >> coord.y = linspace(0,  1, N.y).';                                    %
% >> coord.t = linspace(0, 20, N.t).';                                    %
% >> A       = 0.1;                                                       %
% >> epsi    = 0.25;                                                      %
% >> omega   = 2*pi/10;                                                   %
% >> vel     = doubleGyre(coord, A, epsi, omega);                         %
% >> X       = [reshape(vel.u, [], N.t); reshape(vel.v, [], N.t)];        %
% >>                                                                      %
% >> [Phi, A, props] = podSpaceOnly(X, 'zeromean');                       %
% >>                                                                      %
% >> phi1.u = reshape(Phi(1:N.x*N.y,1), [N.x N.y]);                       %
% >> phi1.v = reshape(Phi(N.x*N.y+1:end,1), [N.x N.y]);                   %
% >> pcolor(coord.x, coord.y, phi1.u.');                                  %
% >> shading interp;                                                      %
% >> axis equal tight;                                                    %
% >> set(gca, 'Layer', 'top');                                            %
% >>                                                                      %
% >> plot(coord.t, A(1,:), 'k');                                          %
% >>                                                                      %
% >> stem(1:20, props.E(1:20)/props.E(1), 'k');                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Phi, A, props] = podSpaceOnly(X, varargin)


%% PARSE INPUTS

% Input defaults.
default.remAvg  = 0;
default.rankTol = 0;

% Input checks.
check.X       = @(x) validateattributes(x,                              ...
                     {'logical', 'numeric'},                            ...
                     {'finite', 'nonempty', 'ndims', 2});
check.remAvg  = @(x) any(validatestring(x, {'zeromean'}));
check.rankTol = @(x) validateattributes(x,                              ...
                     {'logical', 'numeric'},                            ...
                     {'finite', 'nonnegative', 'real'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser , 'X'                         , check.X       );
addOptional ( hParser , 'remAvg'  , default.remAvg  , check.remAvg  );
addParameter( hParser , 'rankTol' , default.rankTol , check.rankTol );
parse(hParser, X, varargin{:});
clear check default;

% Additional verifications.
narginchk(1,4);
nargoutchk(0,3);


%% INITIALIZATIONS

% Determine the number of measurements and realizations.
N.x = size(X, 1);
N.t = size(X, 2);


%% SPACE-ONLY PROPER ORTHOGONAL DECOMPOSITION

% Remove the average.
if hParser.Results.remAvg, X = X - mean(X,2); end

% Solve the eigenvalue problem for the correlation matrix.
[W, E] = eig(X'*X, 'vector');

% Sort the eigenvalues and eigenvectors in descending order.
[E, I] = sort(E, 'descend');
W      = W(:,I.');
clear I;

% Apply the reduced rank approximation.
rk = min(N.x, N.t);
if hParser.Results.rankTol
    E(E < hParser.Results.rankTol*E(1)) = [];

    rk = length(E);
    W  = W(:,1:rk);
end

% Compute the proper orthogonal modes.
Phi = X*W;
Phi = Phi./vecnorm(Phi);
clear W;

% Compute the time dynamics.
A = Phi'*X;
clear X;

% Store decomposition properties.
props.E  = E;
props.H  = abs(-sum((E/sum(E)).*log((E/sum(E))))/log(rk));
props.rk = rk;
clear E rk;


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
% 2023/09/19 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% ~ Add more examples in the preamble.                                    %
% ~ Add other rank reduction options (e.g., based on the total energy).   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
