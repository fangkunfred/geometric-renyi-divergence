% Last Modified: 14 March 2019
% Computes the Rains information of a GADC
% channel, when the input state is sqrt(1-p) |00> + sqrt(p) |11>.
% Written by Khatri et al. in [arXiv:1903.07747]

function [out] = REE_GADC_channel(g,N,p)

Choi_p = [(1-p)*(1-N*g) 0 0 sqrt((1-p)*p*(1-g));
            0 g*N*(1-p) 0 0 ;
            0 0 (1-N)*p*g 0;
            sqrt((1-p)*p*(1-g)) 0 0 p - (1-N)*p*g];

out = REE_qubit_X(Choi_p);

end

