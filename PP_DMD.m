%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                        DYNAMIC MODE DECOMPOSITION                       %
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
% DMD = PP_DMD(VEC, dt);                                                  %
% DMD = PP_DMD(VEC, dt, options);                                         %
% DMD = PP_DMD(VEC, dt, rank, options);                                   %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the dynamic mode decomposition of velocity field %
% data. Several different methods may be selected from, including the     %
% direct method of Rowley et al. (2009), the QR and SVD factorization     %
% methods presented in Schmid (2010), the snapshot method, and the exact  %
% methods presented in Tu et al. (2014).                                  %
%                                                                         %
% References:                                                             %
% Coming Soon ...                                                         %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'DMD'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - Full description coming soon ...                         %
% ----------------------------------------------------------------------- %
% 'dt'         - REAL SCALAR                                              %
%              - Time step between snapshots.                             %
% ----------------------------------------------------------------------- %
% 'options'    - SPECIFIC STRING                                          %
%              - The user may stack options. The mean of the series using %
%                'meandev'. Mode sorting options include 'absFrequency',  %
%                'absRitzFrequency', 'absRitzVal', 'Coherence',           %
%                'Frequency', 'GrowthRate', 'ModeNorm', 'OptimalAmp',     %
%                'OptimalAmpPen', 'RitzFrequency'.  Calculation options   %
%                include 'qr', 'direct', 'exact', 'exactseq', 'snapshots',%
%                and 'svd'.                                               %
% ----------------------------------------------------------------------- %
% 'rank'       - INTEGER SCALAR                                           %
%              - Maximum rank to be used in computing the dynamic mode    %
%                decomposition.                                           %
% ----------------------------------------------------------------------- %
% 'VEC'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - One-dimensional cell array of structs each containing    %
%                information on the spatial mask (C), velocity components %
%                (U and V), and Cartesian grid (X and Y).                 %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Coming Soon ...                                                         %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% field2mat3D                                                             %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_DMD
function [DMD] = PP_DMD(VEC, dt, varargin)

% Determine if the user has requested the options: 1) mean deviation, 2)
% sorting method, 3) DMD method, 4) rank specification.
options = zeros(4,1);
if nargin > 2 && nargin <= 6
    sortopts = { 'absFrequency', 'absRitzFrequency', 'absRitzVal'     , ...
                 'Coherence'   , 'Frequency'       , 'GrowthRate'     , ...
                 'ModeNorm'    , 'OptimalAmp'      , 'OptimalAmpPen'  , ...
                 'RitzFrequency'                                     };
    calcopts = {'qr', 'direct', 'exact', 'exactseq', 'snapshots', 'svd'};
    for k = 1:length(varargin)
        if ischar(varargin{k})
            if strcmpi(varargin{k},'meandev'),     options(1) = k; end
            if max(strcmpi(varargin{k},sortopts)), options(2) = k; end
            if max(strcmpi(varargin{k},calcopts)), options(3) = k; end
        elseif isnumeric(varargin{k}) && isreal(varargin{k})
            if floor(varargin{k}) == varargin{k},  options(4) = k; end
        end
    end
    if ~max(options)
        error('InputError:inOpts',                                      ...
              'Input option(s) not valid.');
    end
else
    if nargin ~= 2, error('InputError:extraOpts',                       ...
                          'Too many input arguments.'); end
end

% Determine the number of grid points in the x and y directions.
[ny, nx] = size(VEC{1}.U);

% Determine the number of time steps (n).
n = length(VEC);

% Reshape the velocity field matrix so that each column represents the
% entire velocity field at a given snapshot.
X = [reshape(field2mat3D(VEC,'U'),[],n);                                ...
     reshape(field2mat3D(VEC,'V'),[],n)];

% If the DMD is to be done with the velocity fields in mean-deviation form,
% subtract off the mean velocity field.
if options(1)
    X = X - repmat(mean(X,2), 1, n);
%     X = X - mean(X,2); % R2017a
end

if options(3) && strcmpi(varargin{options(3)}, 'direct')
    
    % Compute the coefficient vector (s) in the companion matrix (S).
    s = X(:,1:n-1)\X(:,n);
    
    % Construct the companion matrix S.
    S = [zeros(1,n-2) s(1); eye(n-2) s(2:end)];
    
    % Compute the eigenvalues (Mu) of S.
    [~, Mu] = eig(S, 'vector');
    
    % Determine the rank and apply it if necessary.
    if options(4)
        rk = varargin{options(4)};
        Mu(rk+1:end) = 0;
    else
        rk = length(Mu);
    end
    
    % Compute the scaled eigenvectors (Phi) of S.
    Phi = fliplr(vander(Mu));
    if options(4)
        Mu  = Mu(1:rk);
        Phi = Phi(1:rk,:);
    end
    
    % Compute the dynamic modes (Psi).
    Psi = X(:,1:n-1)/Phi;
    
    % Compute the logarithmic mapping (Lambda) of the Ritz values (Mu).
    Lambda = log(Mu)/dt;
    
    % Compute the frequencies of the modes in Hz.
    f = imag(Lambda)/(2*pi);
    
    % Compute the optimal DMD amplitudes.
    A = vecnorm(Psi).';
    
    % Compute the coherence of the modes.
    C = A.*((vecnorm(pinv(Phi)).').^(-1));
    
    % Compute the norms of the dynamic modes.
    N = vecnorm(Psi).';
    
    % Compute a linear amplitude approximation for data reconstruction.
    L = (Psi\X(:,1:n-1)).';
    
    % Sort the modes and computed quantities according to the method
    % specified by the user.
    if options(2)
        switch lower(varargin{options(2)})
            case 'absritzval'
                [~,I] = sort(abs(Mu), 'descend');
            case 'ritzfrequency'
                [~,I] = sort(imag(Mu), 'ascend');
            case 'absritzfrequency'
                [~,I] = sort(abs(imag(Mu)), 'ascend');
            case 'frequency'
                [~,I] = sort(f, 'ascend');
            case 'absfrequency'
                [~,I] = sort(abs(f), 'ascend');
            case 'coherence'
                [~,I] = sort(C, 'descend');
            case 'optimalamp'
                [~,I] = sort(A, 'descend');
            case 'optimalamppen'
                [~,I] = sort(A.*(Mu.^(rk+1)), 'descend');
            case 'growthrate'
                [~,I] = sort(real(Lambda), 'descend');
            case 'modenorm'
                [~,I] = sort(N, 'descend');
        end
        
        A      = A(I);
        C      = C(I);
        f      = f(I);
        L      = L(:,I.');
        Lambda = Lambda(I);
        Mu     = Mu(I);
        N      = N(I);
        Psi    = Psi(:,I.');
    end
    
    % Reshape the dynamic modes to fit the original input dimensions.
    modeU = reshape(Psi(1:nx*ny    ,:), ny, nx, rk);
    modeV = reshape(Psi(nx*ny+1:end,:), ny, nx, rk);
    
    % Store the information of the dynamic mode decomposition in DMD.
    DMD = cell(rk,1);
    for k = 1:rk
        DMD{k} = struct('U'               , modeU(:,:,k),               ...
                        'V'               , modeV(:,:,k),               ...
                        'RitzValues'      , Mu(k)       ,               ...
                        'LogRitzValues'   , Lambda(k)   ,               ...
                        'Frequency'       , f(k)        ,               ...
                        'Coherence'       , C(k)        ,               ...
                        'OptimalAmplitude', A(k)        ,               ...
                        'ModeL2Norm'      , N(k)        ,               ...
                        'LinearApprox'    , L(:,k));
    end
    
elseif options(3) && strcmpi(varargin{options(3)}, 'exact')
    
    % Compute the singular value decomposition of X_1^{n-1}.
    [U, E, V] = svd(X(:,1:n-1),0);
    
    % Determine the rank and apply it if necessary.
    if options(4)
        rk = varargin{options(4)};
        U = U(:,1:rk);
        E = E(1:rk,1:rk);
        V = V(:,1:rk);
    else
        rk = min(size(E));
    end
    
    % Compute the similarity transformation (St) of the companion matrix.
    St = U'*X(:,2:n)*V/E;
    
    % Determine the eigenvalues (Mu) and eigenvectors (Phit) of St.
    [Phit, Mu] = eig(St, 'vector');
    
    % Compute the dynamic modes (Psi).
    Psi = X(:,2:n)*(V/E)*(Phit/diag(Mu));
%     Psi = X(:,2:n)*(V/E)*Phit;
    
    % Compute the logarithmic mapping (Lambda) of the Ritz values (Mu).
    Lambda = log(Mu)/dt;
    
    % Compute the frequencies (f) of the modes in Hz.
    f = imag(Lambda)/(2*pi);
    
    % Compute the coherence of the modes.
    C = (vecnorm(V/E*Phit).^(-1)).';
    
    % Compute the optimal DMD amplitudes.
    A = diag(Mu)\(Psi\X(:,2));
%     A = Psi\X(:,1);
    
    % Compute the norms of the dynamic modes.
    N = vecnorm(Psi).';
    
    % Compute a linear amplitude approximation for data reconstruction.
    L = (Psi\X(:,1:n-1)).';
    
    % Sort the modes and computed quantities according to the method
    % specified by the user.
    if options(2)
        switch lower(varargin{options(2)})
            case 'absritzval'
                [~,I] = sort(abs(Mu), 'descend');
            case 'ritzfrequency'
                [~,I] = sort(imag(Mu), 'ascend');
            case 'absritzfrequency'
                [~,I] = sort(abs(imag(Mu)), 'ascend');
            case 'frequency'
                [~,I] = sort(f, 'ascend');
            case 'absfrequency'
                [~,I] = sort(abs(f), 'ascend');
            case 'coherence'
                [~,I] = sort(C, 'descend');
            case 'optimalamp'
                [~,I] = sort(A, 'descend');
            case 'optimalamppen'
                [~,I] = sort(A.*(Mu.^(rk+1)), 'descend');
            case 'growthrate'
                [~,I] = sort(real(Lambda), 'descend');
            case 'modenorm'
                [~,I] = sort(N, 'descend');
        end
        
        A      = A(I);
        C      = C(I);
        f      = f(I);
        L      = L(:,I.');
        Lambda = Lambda(I);
        Mu     = Mu(I);
        N      = N(I);
        Psi    = Psi(:,I.');
    end
    
    % Reshape the dynamic modes to fit the original input dimensions.
    modeU = reshape(Psi(1:nx*ny    ,:), ny, nx, rk);
    modeV = reshape(Psi(nx*ny+1:end,:), ny, nx, rk);
    
    % Store the information of the dynamic mode decomposition in DMD.
    DMD = cell(rk,1);
    for k = 1:rk
        DMD{k} = struct('U'               , modeU(:,:,k),               ...
                        'V'               , modeV(:,:,k),               ...
                        'RitzValues'      , Mu(k)       ,               ...
                        'LogRitzValues'   , Lambda(k)   ,               ...
                        'Frequency'       , f(k)        ,               ...
                        'Coherence'       , C(k)        ,               ...
                        'OptimalAmplitude', A(k)        ,               ...
                        'ModeL2Norm'      , N(k)        ,               ...
                        'LinearApprox'    , L(:,k));
    end
    
elseif options(3) && strcmpi(varargin{options(3)}, 'exactseq')
    % Exact sequential method of Tu et al. (2014).
    
    % Compute the singular value decomposition of X_1^{n-1}.
    [U, E, V] = svd(X(:,1:n-1),0);
    
    % Determine the rank and apply it if necessary.
    if options(4)
        rk = varargin{options(4)};
        U = U(:, 1:rk);
        E = E(1:rk,1:rk);
        V = V(:, 1:rk);
    else
        rk = min(size(E));
    end
    
    % Compute the similarity transformation (St) of the companion matrix.
    St = U'*X(:,2:n)*V/E;
    
    % Determine the eigenvalues (Mu) and eigenvectors (Phit) of St.
    [Phit, Mu] = eig(St, 'vector');
    
    % Compute the dynamic modes (Psi).
    p = X(:,n) - U*U'*X(:,n);
    if norm(p) < 10^-8
        Psi = U*Phit;
    else
        q = normc(p);
        Psi = U*Phit + q*q'*X(:,2:n)*(V/E)*(Phit/diag(Mu));
    end
    
    % Compute the logarithmic mapping (Lambda) of the Ritz values (Mu).
    Lambda = log(Mu)/dt;
    
    % Compute the frequencies (f) of the modes in Hz.
    f = imag(Lambda)/(2*pi);
    
    % Compute the coherence of the modes.
    C = (vecnorm(V/E*Phit).^(-1)).';
    
    % Compute the optimal DMD amplitudes.
    A = diag(Mu)\(Psi\X(:,2));
    
    % Compute the norms of the dynamics modes.
    N = vecnorm(Psi).';
    
    % Compute a linear amplitude approximation for data reconstruction.
    L = (Psi\X(:,1:n-1)).';
    
    % Sort the modes and computed quantities according to the method
    % specified by the user.
    if options(2)
        switch lower(varargin{options(2)})
            case 'absritzval'
                [~,I] = sort(abs(Mu), 'descend');
            case 'ritzfrequency'
                [~,I] = sort(imag(Mu), 'ascend');
            case 'absritzfrequency'
                [~,I] = sort(abs(imag(Mu)), 'ascend');
            case 'frequency'
                [~,I] = sort(f, 'ascend');
            case 'absfrequency'
                [~,I] = sort(abs(f), 'ascend');
            case 'coherence'
                [~,I] = sort(C, 'descend');
            case 'optimalamp'
                [~,I] = sort(A, 'descend');
            case 'optimalamppen'
                [~,I] = sort(A.*(Mu.^(rk+1)), 'descend');
            case 'growthrate'
                [~,I] = sort(real(Lambda), 'descend');
            case 'modenorm'
                [~,I] = sort(N, 'descend');
        end
        
        A      = A(I);
        C      = C(I);
        f      = f(I);
        L      = L(:,I.');
        Lambda = Lambda(I);
        Mu     = Mu(I);
        N      = N(I);
        Psi    = Psi(:,I.');
    end
    
    % Reshape the dynamic modes to fit the original input dimensions.
    modeU = reshape(Psi(1:nx*ny    ,:), ny, nx, rk);
    modeV = reshape(Psi(nx*ny+1:end,:), ny, nx, rk);
    
    % Store the information of the dynamic mode decomposition in DMD.
    DMD = cell(rk,1);
    for k = 1:rk
        DMD{k} = struct('U'               , modeU(:,:,k),               ...
                        'V'               , modeV(:,:,k),               ...
                        'RitzValues'      , Mu(k)       ,               ...
                        'LogRitzValues'   , Lambda(k)   ,               ...
                        'Frequency'       , f(k)        ,               ...
                        'Coherence'       , C(k)        ,               ...
                        'OptimalAmplitude', A(k)        ,               ...
                        'ModeL2Norm'      , N(k)        ,               ...
                        'LinearApprox'    , L(:,k));
    end
    
elseif options(3) && strcmpi(varargin{options(3)}, 'qr')
    % Compute the QR factorization of X_1^{n-1}.
    [Q, R] = qr(X(:,1:n-1),0);
    
    % Compute the companion matrix S.
    S = (R^(-1))*Q'*X(:,2:n);
    
    % Compute the eigenvalues (Mu) and eigenvectors (Phi) of S.
    [Phi, Mu] = eig(S, 'vector');
    
    % Determine the rank and apply it if necessary.
    if options(4)
        rk  = varargin{options(4)};
        Mu  = Mu(1:rk);
        Phi = Phi(:,1:rk);
    else
        rk = length(Mu);
    end
    
    % Compute the dynamic modes (Psi).
    Psi = X(:,1:n-1)*Phi;
    
    % Compute the logarithmic mapping (Lambda) of the Ritz values (Mu).
    Lambda = log(Mu)/dt;
    
    % Compute the frequencies of the modes in Hz.
    f = imag(Lambda)/(2*pi);
    
    % Compute the coherence of the modes.
    C = vecnorm(Psi).';
    
    % Compute the optimal DMD amplitudes.
    if options(4)
        T = fliplr(vander([Mu;zeros(n-1-rk,1)]));
        T(rk+1:end,:) = [];
        P = diag(C.^(-1))*(Psi'*Psi)*diag(C.^(-1)).*conj(T*T');
        q = conj(diag(T*X(:,1:n-1)'*Psi*diag(C.^(-1))));
        A = P\q;
    else
        T = fliplr(vander(Mu));
        P = diag(C.^(-1))*(Psi'*Psi)*diag(C.^(-1)).*conj(T*T');
        q = conj(diag(T*X(:,1:n-1)'*Psi*diag(C.^(-1))));
        A = P\q;
    end
    
    % Compute the norms of the dynamic modes.
    N = vecnorm(Psi).';
    
    % Compute a linear amplitude approximation for data reconstruction.
    L = (Psi\X(:,1:n-1)).';
    
    % Sort the modes and computed quantities according to the method
    % specified by the user.
    if options(2)
        switch lower(varargin{options(2)})
            case 'absritzval'
                [~,I] = sort(abs(Mu), 'descend');
            case 'ritzfrequency'
                [~,I] = sort(imag(Mu), 'ascend');
            case 'absritzfrequency'
                [~,I] = sort(abs(imag(Mu)), 'ascend');
            case 'frequency'
                [~,I] = sort(f, 'ascend');
            case 'absfrequency'
                [~,I] = sort(abs(f), 'ascend');
            case 'coherence'
                [~,I] = sort(C, 'descend');
            case 'optimalamp'
                [~,I] = sort(A, 'descend');
            case 'optimalamppen'
                [~,I] = sort(A.*(Mu.^(rk+1)), 'descend');
            case 'growthrate'
                [~,I] = sort(real(Lambda), 'descend');
            case 'modenorm'
                [~,I] = sort(N, 'descend');
        end
        
        A      = A(I);
        C      = C(I);
        f      = f(I);
        L      = L(:,I.');
        Lambda = Lambda(I);
        Mu     = Mu(I);
        N      = N(I);
        Psi    = Psi(:,I.');
    end
    
    % Reshape the dynamic modes to fit the original input dimensions.
    modeU = reshape(Psi(1:nx*ny    ,:), ny, nx, rk);
    modeV = reshape(Psi(nx*ny+1:end,:), ny, nx, rk);
    
    % Store the information of the dynamic mode decomposition in DMD.
    DMD = cell(rk,1);
    for k = 1:rk
        DMD{k} = struct('U'               , modeU(:,:,k),               ...
                        'V'               , modeV(:,:,k),               ...
                        'RitzValues'      , Mu(k)       ,               ...
                        'LogRitzValues'   , Lambda(k)   ,               ...
                        'Frequency'       , f(k)        ,               ...
                        'Coherence'       , C(k)        ,               ...
                        'OptimalAmplitude', A(k)        ,               ...
                        'ModeL2Norm'      , N(k)        ,               ...
                        'LinearApprox'    , L(:,k));
    end
    
elseif options(3) && strcmpi(varargin{options(3)}, 'snapshots')
    % Compute the diagonalization of the data.
    [V, E] = eig(X(:,1:n-1)'*X(:,1:n-1));
    
    % Determine the rank and apply it if necessary.
    if options(4)
        rk = varargin{options(4)};
        E = E(1:rk,1:rk);
        V = V(:,1:rk);
    else
        rk = min(size(E));
    end
    
    % Compute the similarity transformation (St) of the companion matrix.
    St = (V/sqrt(E))'*X(:,1:n-1)'*X(:,2:n)*(V/sqrt(E));
    
    % Compute the eigenvalues (Mu) and eigenvectors (Phi) of St.
    [Phit, Mu] = eig(St, 'vector');
    
    % Compute the dynamic modes (Psi).
    Psi = X(:,1:n-1)*(V/sqrt(E))*Phit;
    
    % Compute the logarithmic mapping (Lambda) of the Ritz values (Mu).
    Lambda = log(Mu)/dt;
    
    % Compute the frequencies of the modes in Hz.
    f = imag(Lambda)/(2*pi);
    
    % Compute the coherence of the modes.
    C = (vecnorm(V/sqrt(E)*Phit).^(-1)).';
    
    % Compute the optimal DMD amplitudes.
    if options(4)
        T = fliplr(vander([Mu;zeros(n-1-rk,1)]));
        T(rk+1:end,:) = [];
        A = ((Phit'*Phit).*conj(T*T'))\conj(diag(T*V*sqrt(E)'*Phit));
    else
        T = fliplr(vander(Mu));
        A = ((Phit'*Phit).*conj(T*T'))\conj(diag(T*V*sqrt(E)'*Phit));
    end
    
    % Compute the norms of the dynamic modes.
    N = vecnorm(Psi).';
    
    % Compute a linear amplitude approximation for data reconstruction.
    L = (Psi\X(:,1:n-1)).';
    
    % Sort the modes and computed quantities according to the method
    % specified by the user.
    if options(2)
        switch lower(varargin{options(2)})
            case 'absritzval'
                [~,I] = sort(abs(Mu), 'descend');
            case 'ritzfrequency'
                [~,I] = sort(imag(Mu), 'ascend');
            case 'absritzfrequency'
                [~,I] = sort(abs(imag(Mu)), 'ascend');
            case 'frequency'
                [~,I] = sort(f, 'ascend');
            case 'absfrequency'
                [~,I] = sort(abs(f), 'ascend');
            case 'coherence'
                [~,I] = sort(C, 'descend');
            case 'optimalamp'
                [~,I] = sort(A, 'descend');
            case 'optimalamppen'
                [~,I] = sort(A.*(Mu.^(rk+1)), 'descend');
            case 'growthrate'
                [~,I] = sort(real(Lambda), 'descend');
            case 'modenorm'
                [~,I] = sort(N, 'descend');
        end
        
        A      = A(I);
        C      = C(I);
        f      = f(I);
        L      = L(:,I.');
        Lambda = Lambda(I);
        Mu     = Mu(I);
        N      = N(I);
        Psi    = Psi(:,I.');
    end
    
    % Reshape the dynamic modes to fit the original input dimensions.
    modeU = reshape(Psi(1:nx*ny    ,:), ny, nx, rk);
    modeV = reshape(Psi(nx*ny+1:end,:), ny, nx, rk);
    
    % Store the information of the dynamic mode decomposition in DMD.
    DMD = cell(rk,1);
    for k = 1:rk
        DMD{k} = struct('U'               , modeU(:,:,k),               ...
                        'V'               , modeV(:,:,k),               ...
                        'RitzValues'      , Mu(k)       ,               ...
                        'LogRitzValues'   , Lambda(k)   ,               ...
                        'Frequency'       , f(k)        ,               ...
                        'Coherence'       , C(k)        ,               ...
                        'OptimalAmplitude', A(k)        ,               ...
                        'ModeL2Norm'      , N(k)        ,               ...
                        'LinearApprox'    , L(:,k));
    end
    
elseif ~options(3) || strcmpi(varargin{options(3)}, 'svd')
    % Compute the singular value decomposition of X_1^{n-1}.
    [U, E, V] = svd(X(:,1:n-1),0);
    
    % Determine the rank and apply it if necessary.
    if options(4)
        rk = varargin{options(4)};
        U = U(:, 1:rk);
        E = E(1:rk,1:rk);
        V = V(:, 1:rk);
    else
        rk = min(size(E));
    end
    
    % Compute the similarity transformation (St) of the companion matrix.
    St = U'*X(:,2:n)*V/E;
    
    % Determine the eigenvalues (Mu) and eigenvectors (Phit) of St.
    [Phit, Mu] = eig(St, 'vector');
    
    % Compute the dynamic modes (Psi).
    Psi = U*Phit;
    
    % Compute the logarithmic mapping (Lambda) of the Ritz values (Mu).
    Lambda = log(Mu)/dt;
    
    % Compute the frequencies (f) of the modes in Hz.
    f = imag(Lambda)/(2*pi);
    
    % Compute the coherence of the modes.
    C = (vecnorm(V/E*Phit).^(-1)).';
    
    % Compute the optimal DMD amplitudes.
    if options(4)
        T = fliplr(vander([Mu;zeros(n-1-rk,1)]));
        T(rk+1:end,:) = [];
        A = ((Phit'*Phit).*conj(T*T'))\conj(diag(T*V*E*Phit));
    else
        T = fliplr(vander(Mu));
        A = ((Phit'*Phit).*conj(T*T'))\conj(diag(T*V*E*Phit));
    end
    
    % Compute the norms of the dynamic modes.
    N = vecnorm(Psi).';
    
    % Compute a linear amplitude approximation for data reconstruction.
    L = (Psi\X(:,1:n-1)).';
    
    % Sort the modes and computed quantities according to the method
    % specified by the user.
    if options(2)
        switch lower(varargin{options(2)})
            case 'absritzval'
                [~,I] = sort(abs(Mu), 'descend');
            case 'ritzfrequency'
                [~,I] = sort(imag(Mu), 'ascend');
            case 'absritzfrequency'
                [~,I] = sort(abs(imag(Mu)), 'ascend');
            case 'frequency'
                [~,I] = sort(f, 'ascend');
            case 'absfrequency'
                [~,I] = sort(abs(f), 'ascend');
            case 'coherence'
                [~,I] = sort(C, 'descend');
            case 'optimalamp'
                [~,I] = sort(A, 'descend');
            case 'optimalamppen'
                [~,I] = sort(A.*(Mu.^(rk+1)), 'descend');
            case 'growthrate'
                [~,I] = sort(real(Lambda), 'descend');
            case 'modenorm'
                [~,I] = sort(N, 'descend');
        end
        
        A      = A(I);
        C      = C(I);
        f      = f(I);
        L      = L(:,I.');
        Lambda = Lambda(I);
        Mu     = Mu(I);
        N      = N(I);
        Psi    = Psi(:,I.');
    end
    
    % Reshape the dynamic modes to fit the original input dimensions.
    modeU = reshape(Psi(1:nx*ny    ,:), ny, nx, rk);
    modeV = reshape(Psi(nx*ny+1:end,:), ny, nx, rk);
    
    % Store the information of the dynamic mode decomposition in DMD.
    DMD = cell(rk,1);
    for k = 1:rk
        DMD{k} = struct('U'               , modeU(:,:,k),               ...
                        'V'               , modeV(:,:,k),               ...
                        'RitzValues'      , Mu(k)       ,               ...
                        'LogRitzValues'   , Lambda(k)   ,               ...
                        'Frequency'       , f(k)        ,               ...
                        'Coherence'       , C(k)        ,               ...
                        'OptimalAmplitude', A(k)        ,               ...
                        'ModeL2Norm'      , N(k)        ,               ...
                        'LinearApprox'    , L(:,k));
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
