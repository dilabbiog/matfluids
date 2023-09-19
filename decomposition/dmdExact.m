%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       MODAL DECOMPOSITION TOOLBOX                       %
%                                                                         %
% dmdExact                                                                %
% Dynamic Mode Decomposition                                              %
% Exact dynamic mode decomposition method                                 %
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
% [Phi, A, props] = dmdExact(X);                                          %
% [Phi, A, props] = dmdExact(X, dt);                                      %
% [Phi, A, props] = dmdExact(___, Name, Value);                           %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Compute the exact dynamic mode decomposition [1,2] of a dataset. The    %
% input data is stored in a matrix X whose columns represent distinct     %
% realizations of a measurement (e.g., snapshots over time) and whose     %
% rows represent a series of measurements (e.g., taken at a given time).  %
% The spacing between realizations can be specified (e.g., dt in time or  %
% dx in space). The decomposition can be computed at a reduced rank by    %
% specifying a tolerance (a cutoff threshold based on the fraction of the %
% largest singular value).                                                %
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
% [1] Tu, J. H., Rowley, C. W., Luchtenburg, D. M., Brunton, S. L., &     %
%     Kutz, J. N. (2014). On dynamic mode decomposition: Theory and       %
%     applications. Journal of Computational Dynamics, 1(2), 391-421.     %
% [2] Kutz, J. N., Brunton, S. L., Brunton, B. W., & Proctor, J. L.       %
%     (2016). Dynamic Mode Decomposition: Data-Driven Modeling of Complex %
%     Systems (1st ed.). Society for Industrial and Applied Mathematics   %
%     (SIAM).                                                             %
% [3] Schmid, P. J. (2010). Dynamic mode decomposition of numerical and   %
%     experimental data. Journal of Fluid Mechanics, 656, 5-28.           %
%                                                                         %
% ======================================================================= %
% Input Arguments (Required):                                             %
% ----------------------------------------------------------------------- %
% 'X'            LOGICAL/NUMERIC 2-DIMENSIONAL ARRAY                      %
%              ~ Dataset on which to compute the exact dynamic mode       %
%                decomposition. In the case of fluid flow, the columns of %
%                X often represent sequential snapshots while the rows    %
%                represent the velocity components at various points of   %
%                the flow (within a given snapshot).                      %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'dt'           POSITIVE REAL SCALAR                                     %
%                Default: 1                                               %
%              ~ Spacing between different realizations. In many cases,   %
%                dt is simply the time between sequential snapshots of a  %
%                dataset.                                                 %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% 'pairTol'      NONNEGATIVE REAL SCALAR                                  %
%                Default: 10^(-8)                                         %
%              ~ Tolerance to distinguish between complex conjugate mode  %
%                pairs and distinct modes.                                %
% ----------------------------------------------------------------------- %
% 'rankTol'      NONNEGATIVE REAL SCALAR                                  %
%                Default: 0                                               %
%              ~ Tolerance to reduce the rank of the algorithm at the     %
%                singular value decomposition level. The tolerance is     %
%                taken as a fraction of the largest singular value and is %
%                used as a cutoff (SV < rankTol*max(SV) are ignored).     %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'Phi'          COMPLEX 2-DIMENSIONAL ARRAY                              %
%              ~ Exact dynamic modes of the dataset. The number of rows   %
%                is identical to that of the input data whereas the       %
%                number of columns depends on the number of complex       %
%                conjugate mode pairs as well as on rank reduction.       %
% ----------------------------------------------------------------------- %
% 'A'            COMPLEX 2-DIMENSIONAL ARRAY                              %
%              ~ Time dynamics of the exact dynamic modes of the dataset. %
%                Each row represents the time dynamics of the mode in the %
%                respective column of Phi. The number of columns is       %
%                therefore identical to that of the input data (i.e., the %
%                number of realizations).                                 %
% ----------------------------------------------------------------------- %
% 'props'        STRUCT ARRAY (1 X 1)                                     %
%              ~ Properties of the exact dynamic mode decomposition. The  %
%                structure array contains the following fields:           %
%                 1) amp     : Mode amplitude                             %
%                 2) ampsc   : Scaled mode amplitude                      %
%                 3) coh     : Mode coherence measure [3]                 %
%                 4) eigCT   : Continuous-time eigenvalues                %
%                 5) eigDT   : Discrete-time eigenvalues                  %
%                 6) energy  : Mode energy                                %
%                 7) entropy : Entropy of the energy distribution         %
%                 8) freq    : Mode frequency                             %
%                 9) pair    : Mode has conjugate pair (1) or not (0)     %
%                10) rank    : Rank of the singular value decomposition   %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Phi, A, props] = dmdExact(X, varargin)


%% PARSE INPUTS

% Input defaults.
default.dt      = 1;
default.pairTol = 10^(-8);
default.rankTol = 0;

% Input checks.
check.X       = @(x) validateattributes(x,                              ...
                     {'logical', 'numeric'},                            ...
                     {'finite', 'nonempty', 'ndims', 2});
check.dt      = @(x) validateattributes(x,                              ...
                     {'logical', 'numeric'},                            ...
                     {'finite', 'positive', 'real'});
check.pairTol = @(x) validateattributes(x,                              ...
                     {'logical', 'numeric'},                            ...
                     {'finite', 'nonnegative', 'real'});
check.rankTol = @(x) validateattributes(x,                              ...
                     {'logical', 'numeric'},                            ...
                     {'finite', 'nonnegative', 'real'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'X'                         , check.X       );
addOptional ( hParser, 'dt'      , default.dt      , check.dt      );
addParameter( hParser, 'pairTol' , default.pairTol , check.pairTol );
addParameter( hParser, 'rankTol' , default.rankTol , check.rankTol );
parse(hParser, X, varargin{:});
clear check default;

% Additional verifications.
narginchk(1,6);
nargoutchk(0,3);


%% INITIALIZATIONS

% Determine the number of measurements and realizations.
N.x = size(X, 1);
N.t = size(X, 2);

%% EXACT DYNAMIC MODE DECOMPOSITION

% Compute the singular value decomposition.
[U, S, V] = svd(X(:,1:N.t-1), 'econ');

% Apply the reduced rank approximation.
rk = min(N.x, N.t-1);
if hParser.Results.rankTol
    S = diag(S);
    S(S < hParser.Results.rankTol*S(1)) = [];

    rk = length(S);
    U  = U(:,1:rk);
    S  = diag(S);
    V  = V(:,1:rk);
end

% Compute dynamic modes and discrete-time eigenvalues.
[W, Lambda] = eig(U'*X(:,2:N.t)*V/S, 'vector');
Phi         = X(:,2:N.t)*((V/S)*(W/diag(Lambda)));
C           = (vecnorm(V/S*W).^(-1)).';
clear S U V W;

% Compute dynamic mode properties.
dt    = hParser.Results.dt;
Omega = log(Lambda)/dt;
f     = imag(Omega)/(2*pi);
b     = diag(Lambda)\(Phi\X(:,2));
clear X;


%% ISOLATE UNIQUE FREQUENCIES

% Check for frequency pairs.
idx = 1;
I   = zeros(rk,1);
lst = (1:rk).';
P   = false(rk,1);
tol = hParser.Results.pairTol;
fu  = uniquetol(abs(f), tol, 'DataScale', 1);
for k = 1:length(fu)
    LIA    = ismembertol(abs(f), fu(k), tol, 'DataScale', 1);
    if fu(k) && (sum(LIA) == 2) && (abs(sum(f(LIA))) <= tol)
        I(idx) = max(sign(f(LIA)).*lst(LIA));
        P(idx) = true;
        idx    = idx + 1;
    else
        I(idx:idx+sum(LIA)-1) = lst(LIA);
        idx                   = idx + sum(LIA);
    end
end
I(idx:end) = [];
P(idx:end) = [];
clear k fu LIA lst tol;

% Rearrange the data in order of increasing unique frequency.
f         = f(I);
[f, tmpI] = sort(f, 'ascend');
I         = I(tmpI);
b         = b(I);
C         = C(I);
Lambda    = Lambda(I);
Omega     = Omega(I);
Phi       = Phi(:,I.');
bpen      = abs(b.*(Lambda.^N.t));
clear tmpI;

% Compute the time dynamics.
T                    = fliplr(vander(Lambda));
T(:,length(I)+1:N.t) = T(:,2).^(length(I):N.t-1);
A                    = diag(b)*T;
clear T;

% Compute modal energy and entropy.
E = zeros(length(I),1);
for k = 1:length(I)
    E(k) = (P(k)+1)*dt*trapz(abs(A(k,:)).^2);
end
H = abs(-sum((E/sum(E)).*log((E/sum(E))))/log(rk));
clear dt I k;

% Store decomposition properties.
props.amp     = b;
props.ampsc   = bpen;
props.coh     = C;
props.eigCT   = Omega;
props.eigDT   = Lambda;
props.energy  = E;
props.entropy = H;
props.freq    = f;
props.pair    = P;
props.rank    = rk;
clear b bpen C f H E Lambda Omega P rk;


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
% ~ Add examples in the preamble.                                         %
% ~ Allow for the specification of optional input 'dt' as a 1D array for  %
%   the case of non-sequential data (i.e., non-uniform spacing). Requires %
%   adjustment of continuous time eigenvalues, energy and entropy.        %
% ~ Allow for a selection between a frequency or wavenumber approach for  %
%   continuous eigenvalues.                                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
