% Increasing axes font size

set(0,'defaultaxesfontsize',12)

% Inputting constants

P = [6.4,0.019,4,0.16*2.2,0.72,0.25,6.1,0.3,0.4744];

% Solving using the PGF method to get the expected outbreak probabilities

fun = @(r)PGFmethodebola(r, P);
x0 = [0,0];
x = fsolve(fun,x0);

% Getting parameter values from the vector of parameters
N = P(1);
alpha = P(2);
beta = P(3);
q = P(4);
h = P(5);
lambda_h = P(6);
phi = P(7);
mu = P(8);
gamma = P(9);

% Computing the probability of either recovering from Ebola or dying and having a secure burial
g = gamma*(1-mu)+mu;

% Computing the Reff matrix
Reff = [N*alpha*beta*q+ h*lambda_h, (1-h)*lambda_h;
    h*(N*q*g + (1-g)*(N*q + phi)), (1-h)*(N*q*g + (1-g)*(N*q + phi))];

K = Reff;

% Initialising constants
beta_hh = K(1,1)/6.5;
beta_hc = K(1,2)/6.5;
beta_ch = K(2,1)/15;
beta_cc = K(2,2)/15;


nsims = 1000000;
% Creating table for the maximum number of infections
maxIvalues = zeros(nsims,1);

for j = 1:nsims

    % Initialising sizes of each population group

    I_h = 1;
    I_c = 0;
    
    I = I_c + I_h;
    R = 0;
    totalI = I+R;

    % Gillespie algorithm loop
    while I>0 && I+R<501
        lambda = (beta_hh+beta_hc)*I_h+(beta_ch+beta_cc)*I_c+(1/6.5)*I_h+(1/15)*I_c;

         r = rand();
         if (r < ((beta_hh*I_h+beta_ch*I_c)/(lambda)))
             I_h = I_h+1;
         elseif (r < ((beta_hh+beta_hc)*I_h+(beta_ch+beta_cc)*I_c)/lambda)
             I_c = I_c+1;
         elseif (r < (lambda - (1/15)*I_c)/lambda)
             I_h = I_h-1;
             R = R+1;
         else
             I_c = I_c-1;
             R = R+1;
         end
         I = I_c + I_h;
    end   
    maxIvalues(j) = I+R;
    
end

% Generating array of proportion of tests where different thresholds are
% exceeded

B = zeros(501,2);
B(:, 1) = 0:500;
    
for i = 1:501
     B(i,2) = (sum(maxIvalues >= B(i,1))/nsims);
end

% Plotting lines

figure(2)
plot(B(:,1), B(:,2), 'Color', 'blue');
hold on
yline(max(0, 1-x(1)), 'Color', 'blue', 'LineStyle','--');

% Starting again with 1 infected community patient initially

maxIvalues = zeros(nsims,1);

for j = 1:nsims

    % Initialising sizes of each population group

    I_h = 0;
    I_c = 1;
    
    I = I_c + I_h;
    R = 0;
    totalI = I+R;

    % Gillespie algorithm loop
    while I>0 && I+R<501
        lambda = (beta_hh+beta_hc)*I_h+(beta_ch+beta_cc)*I_c+(1/6.5)*I_h+(1/15)*I_c;

         r = rand();
         if (r < ((beta_hh*I_h+beta_ch*I_c)/(lambda)))
             I_h = I_h+1;
         elseif (r < ((beta_hh+beta_hc)*I_h+(beta_ch+beta_cc)*I_c)/lambda)
             I_c = I_c+1;
         elseif (r < (lambda - (1/15)*I_c)/lambda)
             I_h = I_h-1;
             R = R+1;
         else
             I_c = I_c-1;
             R = R+1;
         end
         I = I_c + I_h;
    end   
    maxIvalues(j) = I+R;
    
end

% Generating array of proportion of tests where different thresholds are
% exceeded

B = zeros(501,2);
B(:, 1) = 0:500;
    
for i = 1:501
     B(i,2) = (sum(maxIvalues >= B(i,1))/nsims);
end

% Plotting lines

figure(2)
% Plotting lines for the proportion of outbreaks where thresholds for the maximum number of people infected simultaneously are exceeded
plot(B(:,1), B(:,2), 'Color', 'r');
hold on
% Plotting the expected probabilities of major outbreak
yline(max(0, 1-x(2)), 'Color', 'r', 'LineStyle','--');
% Adding labels and legend
xlabel('Total Number of people infected (M)')
ylabel('Proportion of simulations in which \newline threshold M was reached')
legend('Simulation starting from I_h = 1, I_c = 0', 'Expected PMO starting from I_c = 1, I_h = 0', 'Simulation starting from I_h = 0, I_c = 1', 'Expected PMO starting from I_c = 0, I_h = 1')
