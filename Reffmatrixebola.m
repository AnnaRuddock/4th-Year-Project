P = [6.4;0.019;4;0.16*2.2;0.72;0.25;6.1;0.3;0.4744];

N = P(1);
alpha = P(2);
beta = P(3);
q = P(4);
h = P(5);
lambda_h = P(6);
phi = P(7);
mu = P(8);
gamma = P(9);

g = (1-mu)*gamma+mu;

Reff = [N*alpha*beta*q+ h*lambda_h, (1-h)*lambda_h;
    h*(N*q*g + (1-g)*(N*q + phi)), (1-h)*(N*q*g + (1-g)*(N*q + phi))]
