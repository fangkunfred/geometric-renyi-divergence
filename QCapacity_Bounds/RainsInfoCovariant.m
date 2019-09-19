% compute the Rains information for covariant channels
% input JN: the Choi matrix of the channel
%        d: dimension vector in the order of [da db]
% Writte by Kun Fang

function r = RainsInfoCovariant(JN,d)
da = d(1);
db = d(2);
JN = (JN+JN')/2; % avoid numerical problems from non-Hermitian JN

cvx_begin sdp quiet
    variable sig(da*db,da*db) hermitian
    variable X(da*db,da*db) hermitian
    variable Y(da*db,da*db) hermitian
    r = quantum_rel_entr(JN/da,sig)/log(2);
    minimize r
    subject to
        sig >= 0; X >= 0, Y >= 0; 
        PartialTranspose(sig,2,[da db]) == X - Y; trace(X + Y) <= 1;
cvx_end
end