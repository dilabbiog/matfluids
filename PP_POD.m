%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                     PROPER ORTHOGONAL DECOMPOSITION                     %
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
% SYNTAX                                                                  %
%                                                                         %
% [POD,H] = PP_POD(VEC);                                                  %
% [POD,H] = PP_POD(VEC, options);                                         %
% [POD,H] = PP_POD(VEC, rank, options);                                   %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the proper orthogonal decomposition (POD) of     %
% velocity field data using the algorithm described by L. Sirovich [1].   %
% The algorithm is also well-described by K. Meyer [2]. The user has the  %
% option of performing the POD in mean deviation form where the mean flow %
% is subtracted off at each time step. This function also computes the    %
% global entropy of the flow as defined by N. Aubry [3].                  %
%                                                                         %
% References:                                                             %
% [1] Sirovich, L. (1987). Turbulence and the dynamics of coherent        %
%     structures, Part I: Coherent structures. Quarterly of Applied       %
%     Mathematics, 45(3), 561-571.                                        %
% [2] Meyer, K., Pedersen, J. M. & Özkan, O. (2007). A turbulent jet in   %
%     crossflow analysed with proper orthogonal decomposition. Journal of %
%     Fluid Mechanics, 583, 199-227.                                      %
% [3] Aubry, N. (1991). On the hidden beauty of proper orthogonal         %
%     decomposition. Theoretical and Computational Fluid Dynamics, 2,     %
%     339-352.                                                            %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'H'          - REAL SCALAR                                              %
%              - Shannon entropy of the energy distribution from the      %
%                proper orthogonal decomposition. The entropy is a        %
%                measure of energy spreading between the POD modes and    %
%                approaches zero as all the energy concentrates in the    %
%                first mode and conversely approaches unity as the energy %
%                becomes equally distributed among all the modes.         %
% ----------------------------------------------------------------------- %
% 'POD'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - Proper orthogonal decomposition data containing the mode %
%                velocity field components 'U' and 'V', amplitudes 'a',   %
%                and energy characteristics 'E' (with subfields 'Energy', %
%               'Frac', 'CumSum' and 'CumFrac').                          %
% ----------------------------------------------------------------------- %
% 'VEC'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - Velocity field data containing the mask 'C', Cartesian   %
%                coordinates 'X', 'Y' and 'Z', and velocity field         %
%                components 'U', 'V' and 'W'.                             %
% ----------------------------------------------------------------------- %
% Options:                                                                %
% ----------------------------------------------------------------------- %
% 'meandev'    - SPECIFIC STRING                                          %
%              - Option to apply the proper orthogonal decomposition in   %
%                mean deviation form using 'meandev'.                     %
%              - Default: Not Active                                      %
% ----------------------------------------------------------------------- %
% 'rank'       - POSITIVE, NON-ZERO INTEGER                               %
%              - Rank of the proper orthogonal decomposition.             %
%              - Default: length(VEC)                                     %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Perform a proper orthogonal decomposition on a time-dependent double    %
% gyre on the domain (x,y) = [0,2]x[0,1] with a constant grid spacing of  %
% 0.01 over the time interval [0,20] with time-step size 0.1. Use A = 0.1,%
% epsilon = 0.25, and omega = 2*pi/10. Cycle through the modes to observe %
% their structures. Plot the cumulative fraction of total energy vs. mode %
% number. Plot a bar chart of the energy fractions as well.               %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> [POD,H] = PP_POD(VEC);                                               %
% >> for k = 1:length(t)                                                  %
% quiver(VEC{k}.X(1:4:end,1:4:end), VEC{k}.Y(1:4:end,1:4:end), ...        %
%        POD{k}.U(1:4:end,1:4:end), POD{k}.V(1:4:end,1:4:end), 'k');      %
% axis([0 2 0 1]);                                                        %
% pause(0.25);                                                            %
% end                                                                     %
% >> clear k;                                                             %
% >> plot(1:length(t), field2vec(field2cell(POD,'E'),'CumFrac'), 'k');    %
% >> axis([1 length(t) 0 1]);                                             %
% >> bar(1:length(t), field2vec(field2cell(POD,'E'),'Frac'));             %
% >> axis([0.5 length(t) 0 1]);                                           %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% field2mat3D                                                             %
% isint                                                                   %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_POD
function [POD, H] = PP_POD(VEC, varargin)

%% Parse the inputs.

% Determine if the user has requested to compute POD in mean deviation form
% (option 1) and/or if a rank has been specified (option 2).
options = zeros(2,1);
if nargin > 1 && nargin <= 3
    for k = 1:length(varargin)
        if ischar(varargin{k})
            if strcmpi(varargin{k},'meandev')
                options(1) = k;
            else
                error('InputError:inCharOpt',                                      ...
                     ['The specified character or string input option ' ...
                      'is not supported.']);
            end
        elseif isnumeric(varargin{k})
            if isint(varargin{k}) && length(varargin{k}) == 1           ...
                                  && varargin{k} > 0
                varargin{k} = floor(abs(varargin{k}));
                options(2) = k;
            else
                error('InputError:inIntOpt',                                      ...
                     ['The rank must be a single integer value '        ...
                      'greater than zero.']);
            end
        end
    end
    if ~max(options)
        error('InputError:inOpts',                                      ...
              'Input option(s) not valid.');
    end
else
    if nargin ~= 1, error('InputError:extraOpts',                       ...
                          'Too many input arguments.'); end
end

%% Prepare the data for POD.

% Determine the number of data points in x (nx) and y (ny).
[ny, nx] = size(VEC{1}.X);

% Determine the number of time steps (n).
n = length(VEC);

% Reshape the velocity field matrix so that each column represents the
% entire velocity field at a given snapshot.
X = [reshape(field2mat3D(VEC,'U'),[],n); ...
     reshape(field2mat3D(VEC,'V'),[],n)];

%% Perform the mean subtraction.
 
% If the POD is to be done with the velocity fields in mean-deviation form,
% subtract off the mean velocity field.
if options(1), X = X - repmat(mean(X,2), 1, n); end
%              X = X - mean(X,2); % R2017a

%% Compute the POD.

% Compute the covariance matrix (C).
C = X'*X;

% Determine the eigenvalues (S) and eigenvectors (V) of the covariance
% matrix.
[V,S] = eig(C, 'vector');

% Order the eigenvalues (and corresponding eigenvectors) in descending
% order.
[S,I] = sort(S, 'descend');
V = V(:, I.');

% Apply the rank.
if options(2)
    rk = varargin{options(2)};
    S = S(1:rk);
    V = V(:, 1:rk);
else
    rk = length(S);
end

% Compute the normalized proper orthogonal modes.
M = normc(X*V);

% Compute the proper orthogonal mode coefficients (or weights).
a = M'*X;

%% Compute additional POD parameters.

% Compute the fraction of total energy, the cumulative energy, and the
% cumulative fraction of total energy for all time instants.
Efrac  = S/sum(S);
Ecsum  = cumsum(S);
Ecfrac = Ecsum/sum(S);

% Compute the entropy (a measure of energy spread) between the modes.
H = abs(-sum(Efrac.*log(Efrac))/log(rk));

%% Prepare the POD data for output.

% Reshape the proper orthogonal modes to fit the original input dimensions.
modeU = reshape(M(1:nx*ny    ,:), ny, nx, rk);
modeV = reshape(M(nx*ny+1:end,:), ny, nx, rk);

% Initialize a cell array of structs (POD) to hold the proper orthogonal
% modes of the velocity components (U and V), the coefficients (a), and the
% energy (E).
POD = cell(rk,1);
for k = 1:rk
    POD{k} = struct('U', modeU(:,:,k), 'V', modeV(:,:,k),               ...
                    'a', a(k,:).'    , 'E', 0);
    POD{k}.E = struct('Energy', S(k)    , 'Frac'   , Efrac(k),          ...
                      'CumSum', Ecsum(k), 'CumFrac', Ecfrac(k));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% SUPPRESS MESSAGES %%%%%%%%%%%%%%%%%%%%%%%%% %%

%#ok<*N/A>
% Line(s) N/A
% Message(s)
% * N/A.
% Reason(s)
% * N/A.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Line(s) 178, 182, 186-187                                               %
% * Note that although the eigendecomposition method was used, one could  %
%   also use the singular value decomposition method, which would involve %
%   replacing these lines with the two lines:                             %
%   [~, S, V] = svd(X,0);                                                 %
%   S = diag(S).^2;                                                       %
% * Although the method seems much simpler to implement, it will execute  %
%   significantly slower than the eigendecomposition method used here. It %
%   also has significant memory limitations.                              %
%                                                                         %
% Line(s) 199, 202                                                        %
% * It should be noted that the POD modes, being eigenvectors, are bi-    %
%   directional vectors. We can therefore multiply any POD mode by -1 to  %
%   change its direction provided of course its time-varying amplitude is %
%   is also multiplied by -1.                                             %
%                                                                         %
% Line(s) 202                                                             %
% * Note that the POD mode matrix M is a unitary matrix and therefore we  %
%   could have written a = M\X as well, although more computationally     %
%   expensive for large M.                                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
