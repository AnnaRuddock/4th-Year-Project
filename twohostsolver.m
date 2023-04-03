% Importing the PGF method function
fun = @PGFmethod;
x0 = [0,0];
% Changing options so no display shows in the command window when functions are solved
options = optimoptions('fsolve','Display','none');
% Using fsolve to find a vector that satisfied the equations in PGFmethod
x = fsolve(fun,x0, options)

% Importing the first step method function
fun = @firststepmethod;
y0 = [0,0];
% Using fsolve to find a vector that satisfied the equations in firststepmethod
y = fsolve(fun,y0, options)
