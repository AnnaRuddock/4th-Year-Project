para = struct('N', 6.4, 'alpha', 0.019, ...
    'beta', 4, 'q', 0.16*2.2, 'h', 0.72, 'lambda_h', 0.25, 'phi', 6.1, ...
    'mu', 0.3, 'gamma_1', 0.68, 'gamma_2', 1.9, 't', 8);

gamma = para.gamma_1*(1-1/((para.t-7)*para.gamma_2));

g = gamma+para.mu;

Reff = [para.N*para.alpha*para.beta*para.q+ para.h*para.lambda_h, (1-para.h)*para.lambda_h;
    para.h*(para.N*para.q*g + (1-g)*(para.N*para.q + para.phi)), (1-para.h)*(para.N*para.q*g + (1-g)*(para.N*para.q + para.phi))]