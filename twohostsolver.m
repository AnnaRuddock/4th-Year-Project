fun = @PGFmethod;
x0 = [0,0];
options = optimoptions('fsolve','Display','none');
x = fsolve(fun,x0, options)

fun = @firststepmethod;
y0 = [0,0];
y = fsolve(fun,y0, options)