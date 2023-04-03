% Inputting vector of parameter values
P = [6.4,0.019,4,0.16*2.2,0.72,0.25,6.1,0.3,0.4744];

% Solving using the PGF method to get the expected outbreak probabilities

fun = @(r)PGFmethodebola(r, P);
x0 = [0,0];
options = optimoptions('fsolve','Display','none');
x = fsolve(fun,x0, options);
x = 1.-x

% Now solving with the first step analysis method to estimate the same probabilities
fun = @(r)firststepmethodebola(r,P);
y0 = [0 0];
y = fsolve(fun, y0, options);
y = 1.-y
