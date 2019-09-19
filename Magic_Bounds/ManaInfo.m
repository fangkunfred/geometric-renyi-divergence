% compute the mana of a quantum channel
% input channel Choi matrix JN and input output dimension d = [da db]
% Writte by Kun Fang
function t = ManaInfo(JN,d)

X = GenPauli(1,0,d);
Z = GenPauli(0,1,d);
w = exp(2*pi*1j/d);
tau = exp((d+1)*pi*1j/d);

A0 = zeros(d);
for a1 = 1:d
for a2 = 1:d
T(:,:,a1,a2) = tau^(-a1*a2)*Z^(a1)*X^(a2);
A0 = A0 + T(:,:,a1,a2);
end
end
A0 = A0/d;

sum = 0;
x = zeros(d);
for a1 = 1:d
for a2 = 1:d
A(:,:,a1,a2) = T(:,:,a1,a2)*A0*T(:,:,a1,a2)';
x(a1,a2) = WignerTraceNorm(ApplyMap(A(:,:,a1,a2),JN),d);
end
end
t = max(max(x));
t = log2(t);
end