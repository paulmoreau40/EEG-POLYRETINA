function P = prewhiteningMatrix(C, noiselambda, Tol_svd)
% From ft_inverse_mne (around lines 210)
% P corresponds to C^(-1/2)

% Optional parameters
% Default values
if nargin < 3
    % Tolerance for singular value decomposition
    Tol_svd = 1e-12;
end
if nargin < 2
    %Regularization parameter
    noiselambda = 0;
end

% Scale C so that tr(C) = tr(Id)
% (for comparison with the main term)
%C = C.*(size(C,1)/trace(C));
% Compute singular value decomposition
C_reg = C + eye(size(C))*noiselambda;
[U,S,~] = svd(C_reg);
diagS   = diag(S);
sel     = find(diagS>Tol_svd.*diagS(1));
P       = diag(1./sqrt(diag(S(sel,sel))))*U(:,sel)';

% Check that the square root inversion happened correctly
C_reconstructed = eye(size(C_reg))/(P'*P);
Tol_reconstruction = 1000*Tol_svd;
if all((C_reg - C_reconstructed)<Tol_reconstruction, 'all')
    fprintf('Pseudo square root inversion successfull with %.0e regularization (%.0e tolerance)\n',...
        noiselambda, Tol_reconstruction)
else
    warning('Pseudo square root inversion failed with %.0e regularization (%.0e tolerance)',...
        noiselambda, Tol_reconstruction)
end
end