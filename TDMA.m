%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                       TRIDIAGONAL MATRIX ALGORITHM                      %
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
% X = TDMA(a, b, c, F);                                                   %
% X = TDMA(a, b, c, F, 'sparse');                                         %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function applies the tridiagonal matrix algorithm (TDMA), also     %
% known as the Thomas algorithm (attributed to Llewellyn Hilleth Thomas   %
% from [1]), to solve a tridiagonal system of equations.                  %
%                                                                         %
% References:                                                             %
% [1] Thomas, L. H. (1949). Elliptic problems in linear difference        %
%     equations over a network (Tech. Rep.). New York, NY: Watson         %
%     Scientific Computing Laboratory, Columbia University.               %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'a'          - REAL VECTOR                                              %
%              - One-dimensional real scalar array representing the lower %
%                diagonal vector of a tridiagonal matrix having n - 1     %
%                entries for a matrix of size n.                          %
% ----------------------------------------------------------------------- %
% 'b'          - REAL VECTOR                                              %
%              - One-dimensional real scalar array representing the       %
%                diagonal vector of a tridiagonal matrix having n entries %
%                for a matrix of size n.                                  %
% ----------------------------------------------------------------------- %
% 'c'          - REAL VECTOR                                              %
%              - One-dimensional real scalar array representing the upper %
%                diagonal vector of a tridiagonal matrix having n - 1     %
%                entries for a matrix of size n.                          %
% ----------------------------------------------------------------------- %
% 'F'          - REAL VECTOR                                              %
%              - One-dimensional real scalar array representing the       %
%                resultant vector for the system AX = F, having n entries %
%                for a matrix A of size n.                                %
% ----------------------------------------------------------------------- %
% 'x'          - REAL VECTOR                                              %
%              - One-dimensional real scalar array representing the       %
%                soultion vector for the system AX = F, having n entries  %
%                for a matrix A of size n.                                %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Solve the system of n equations AX = F given some random lower, mid,    %
% and upper diagonal matrices 'a', 'b' and 'c' as well as some random     %
% resultant vector 'F'.                                                   %
%                                                                         %
% >> n = 50;                                                              %
% >> a = rand(n - 1, 1);                                                  %
% >> b = rand(n, 1);                                                      %
% >> c = rand(n - 1, 1);                                                  %
% >> F = rand(n, 1);                                                      %
% >> X_TDMA = TDMA(a, b, c, F);                                           %
% >> X_SPRS = TDMA(a, b, c, F, 'sparse');                                 %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% N/A                                                                     %
%                                                                         %
% Called in:                                                              %
% FD_COMP2                                                                %
% FD_COMP3                                                                %
% FD_COMP4                                                                %
% FD_EXO2                                                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TDMA
function [X] = TDMA(a, b, c, F, varargin)

% Determine the number of equations (n) in the system from the length of
% the resultant vector (F).
n = length(F);

if nargin == 4
    % Use the tridiagonal matrix algorithm (TDMA).
    
    % Initialize the vectors to be used in the algorithm, including the
    % solution vector (X).
    beta = zeros(n,1);
    gamma = zeros(n-1,1);
    g = zeros(n,1);
    X = zeros(n,1);
    
    % Perform the LU decomposition.
    beta(1) = b(1);
    g(1) = F(1);
    for k = 2:n
        gamma(k-1) = a(k-1)/beta(k-1);
        beta(k) = b(k) - gamma(k-1)*c(k-1);
        g(k) = F(k) - gamma(k-1)*g(k-1);
    end
    
    % Perform the back-substitution.
    X(n) = g(n)/beta(n);
    for k = n-1:-1:1
        X(k) = (g(k) - c(k)*X(k+1))/beta(k);
    end
    
elseif nargin == 5 && strcmpi(varargin{1}, 'sparse')
    % Use MATLAB's sparse matrix methods.
    
    % Create a sparse tridiagonal matrix (A) with lower diagonal 'a',
    % diagonal 'b', and upper diagonal 'c'.
    A = spdiags([[a; 0] b [0; c]], -1:1, n, n);
    
    % Compute the solution vector using MATLAB's backslash.
    X = A\F;
    
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
% Line(s) 120                                                             %
% * Alternatively, strcmp() could be used (which may improve the speed of %
%   the function if it is within a loop), however keep in mind that the   %
%   comparison will be case sensitive as opposed to strcmpi() which is    %
%   case insensitive. In this case, the call occurs but once and cannot   %
%   be expected to improve computation time significantly.                %
%                                                                         %
% Line(s) 125                                                             %
% * Note the zeros appended to the vectors 'a' and 'c'. The spdiags()     %
%   function requires the diagonal vector inputs to have the same length, %
%   however for lower diagonal number 'p', the function will store the    %
%   first 'n - p' entries and for upper diagonal number 'q', the function %
%   will store the last 'n - q' entries.                                  %
%                                                                         %
% Line(s) 128                                                             %
% * MATLAB's backslash is very efficient as it takes into account several %
%   algorithms for different matrix types. However, in this case, the     %
%   computation time required when using spdiags() combined with MATLAB's %
%   backslash is longer for larger 'n' than when using the tridiagonal    %
%   matrix algorithm.                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    