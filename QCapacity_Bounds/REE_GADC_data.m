% get data and plot for Rains information of a GADC
% Written by Khatri et al. in [arXiv:1903.07747]
clc;clear;

% N_vec = [0,0.3,0.5];
N_vec=[0 0.3 0.45];

g_vec = [0:0.01:1];

tic
for i = 1:length(N_vec)
    i
    N = N_vec(i);
    
    for j = 1:length(g_vec)
        
        j
            
        g = g_vec(j);
        
        [tmp, f] = fminsearch(@(p) -REE_GADC_channel(g,N,p), [1/2]);
        
        Rains_inf(i,j) = -f;
        
        disp(['N value = ' num2str(N) ', g value = ' num2str(g) ', Rains info = ' num2str(-f)]);
        
    end
    
end
toc


