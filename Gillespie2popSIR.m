% clearing all stored values and the command window

clc
clear
close all

% Solving using the PGF method to get the expected outbreak probabilities

fun = @PGFmethod;
x0 = [0,0];
x = fsolve(fun,x0);

% Inputting the NGM matrix

K = [1.41, 0.34; 0.35, 0.87];

% propc = proportion of population that are children
propc = 0.32;

% Initialising constants
u = 1;
N = 1000;
N_c = propc*N;
N_a = N-N_c;
Beta = K*u;
beta_cc = Beta(1,1)/N_c;
beta_ac = Beta(1,2)/N_c;
beta_ca = Beta(2,1)/N_a;
beta_aa = Beta(2,2)/N_a;

nsims = 10000;
% Creating table for the maximum number of infections
maxIvalues = zeros(nsims,1);
    
for j = 1:nsims

   %  if mod(j,100000) == 0
    %        j
     %end

     % Initialising sizes of each population group

     S_c = N_c - 1;
     I_c = 1;
     R_c = 0;
     S_a = N_a;
     I_a = 0;
     R_a = 0;

     I = I_c + I_a;
     maxI = I;

     % Gillespie algorithm loop
     while I>0
        lambda = beta_cc*S_c*I_c + beta_ac*S_c*I_a + beta_ca*S_a*I_c + beta_aa*S_a*I_a + u*(I_c+I_a);

         r = rand();
         if (r < ((beta_cc*S_c*I_c + beta_ac*S_c*I_a)/(lambda)))
             S_c = S_c-1;
             I_c = I_c+1;
         elseif (r < (lambda - u*(I_c + I_a))/lambda)
             S_a = S_a-1;
             I_a = I_a+1;
         elseif (r < (lambda - u*I_a)/lambda)
             I_c = I_c-1;
             R_c = R_c+1;
         else
             I_a = I_a-1;
	         R_a = R_a+1;
         end
         I = I_c + I_a;
         if (I > maxI)
	     maxI = I;
         end
     end
        
     maxIvalues(j) = maxI;
    
end
 
% Generating array of proportion of tests where different thresholds are
% exceeded

B = zeros(151,2);
B(:, 1) = 0:150;
    
for i = 1:151
     B(i,2) = (sum(maxIvalues >= B(i,1))/nsims);
end

% Plotting lines

figure(1)
plot(B(:,1), B(:,2), 'Color', 'blue');
hold on
yline(max(0, 1-x(1)), 'Color', 'blue', 'LineStyle','--');

% starting again with 1 infected adult initially

maxIvalues = zeros(nsims,1);
    
for j = 1:nsims

    % if mod(j,100000) == 0
    %        j
    % end

     S_c = N_c;
     I_c = 0;
     R_c = 0;
     S_a = N_a - 1;
     I_a = 1;
     R_a = 0;

     I = I_c + I_a;
     maxI = I;
     while I>0
        lambda = beta_cc*S_c*I_c + beta_ac*S_c*I_a + beta_ca*S_a*I_c + beta_aa*S_a*I_a + u*(I_c+I_a);

         r = rand();
         if (r < ((beta_cc*S_c*I_c + beta_ac*S_c*I_a)/(lambda)))
             S_c = S_c-1;
             I_c = I_c+1;
         elseif (r < (lambda - u*(I_c + I_a))/lambda)
             S_a = S_a-1;
             I_a = I_a+1;
         elseif (r < (lambda - u*I_a)/lambda)
             I_c = I_c-1;
             R_c = R_c+1;
         else
             I_a = I_a-1;
	         R_a = R_a+1;
         end
         I = I_c + I_a;
         if (I > maxI)
	     maxI = I;
         end
     end
        
     maxIvalues(j) = maxI;
    
end
    
B = zeros(151,2);
B(:, 1) = 0:150;
    
for i = 1:151
     B(i,2) = (sum(maxIvalues >= B(i,1))/nsims);
end

plot(B(:,1), B(:,2), 'Color', 'red');
yline(max(0, 1-x(2)), 'Color', 'red', 'LineStyle','--');

% Adding labels to axes

xlabel('Number of individuals infected simultaneously (M)')
ylabel('Proportion of simulations in which threshold M was reached')

% Adding legend labels

legend({'Start with 1 infected child', '', 'Start with 1 infected adult', ''}, 'Location','northeast')


