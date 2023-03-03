function F = PGFmethodebola(r, P)

N = P(1);
alpha = P(2);
beta = P(3);
q = P(4);
h = P(5);
lambda_h = P(6);
phi = P(7);
mu = P(8);
gamma = P(9);

g = gamma*(1-mu)+mu;

Reff = [N*alpha*beta*q+ h*lambda_h, (1-h)*lambda_h;
    h*(N*q*g + (1-g)*(N*q + phi)), (1-h)*(N*q*g + (1-g)*(N*q + phi))];

K = Reff;

F(1) = r(1) - 1/(K(1,1)*(1-r(1)) + K(1,2)*(1-r(2))+1);
F(2) = r(2) - 1/(K(2,1)*(1-r(1)) + K(2,2)*(1-r(2))+1);

end