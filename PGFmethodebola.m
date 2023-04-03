function F = PGFmethodebola(r, P)

% Getting the parameters from the inputted vector
N = P(1);
alpha = P(2);
beta = P(3);
q = P(4);
h = P(5);
lambda_h = P(6);
phi = P(7);
mu = P(8);
gamma = P(9);

% Computing the probability of either recovering or dying and having a secure burial
g = gamma*(1-mu)+mu;

% Computing the Reff matrix
Reff = [N*alpha*beta*q+ h*lambda_h, (1-h)*lambda_h;
    h*(N*q*g + (1-g)*(N*q + phi)), (1-h)*(N*q*g + (1-g)*(N*q + phi))];

K = Reff;

%Writing the equations that r(1), r(2) (the probabilities of no major outbreak) must solve
F(1) = r(1) - 1/(K(1,1)*(1-r(1)) + K(1,2)*(1-r(2))+1);
F(2) = r(2) - 1/(K(2,1)*(1-r(1)) + K(2,2)*(1-r(2))+1);

end
