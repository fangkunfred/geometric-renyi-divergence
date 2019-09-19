% Last Modified: 14 March 2019
% Calculates the Rains information of a bipartite qubit
% state.
% Written by Khatri et al. in [arXiv:1903.07747]
function [out] = REE_qubit_X(rho_AB,varargin)

cvx_begin sdp

%cvx_solver sedumi

cvx_expert true

%if ~isempty(varargin)
	cvx_quiet true
%end

variable sigma_AB(4,4) hermitian
variables a b c d  
variable xi complex

RE = quantum_rel_entr(rho_AB,sigma_AB)/log(2);

minimize RE

sigma_AB == 1/2 * [a 0 0 xi; 0 b 0 0; 0 0 c 0; conj(xi) 0 0 d];

sigma_AB >= 0;

trace(sigma_AB) == 1;

PartialTranspose(sigma_AB,2,[2,2]) >= 0;

cvx_end

out=full(RE);

opt_state=full(sigma_AB);

end

