clc
clear
close all

% Inputting a range of values for R_0
Rvalues = [0.8, 1.01, 1.1, 1.3, 1.6, 2];
colour = ["r" "g" "b" "c" "m" "k"];


figure(1)
for k = 1:6

    % Setting gamma=1 and hence beta=R_0 (since it is the value of R_0 that matters rather than what beta and gamma are specifically)
    gamma = 1;
    N = 1000;
    beta = (Rvalues(k)*u)/N; %Scaling beta since we assumed populations add up to 1 rather than N
    nsims = 1000000;
    maxIvalues = zeros(nsims,1);
    
    for j = 1:nsims

	% initialising with 1 infected person initially and everyone else susceptible
        S = 999;
        I = 1;
        maxI = I;
        R = 0;
        while I>0
	    % Computing total rate of events lambda
            lambda = ((beta*S*I))+u*I;

	    % Generating a random number r, choosing an event to occur and implementing it
            r = rand();
            if (r < (((beta*S*I))/(((beta*S*I))+u*I)))
                S = S-1;
                I = I+1;
            else
                I = I-1;
	            R = R+1;
            end
            if (I > maxI)
	        maxI = I;
            end
        end
        
	% Storing the maximum number of people infected simultaneously
        maxIvalues(j) = maxI;
    
    end
    
    B = zeros(251,2);
    B(:, 1) = 0:250;
    
    % Computing the proportion of simulations in which thresholds for the maximum number of people infected simultaneously are reached
    for i = 1:251
        B(i,2) = (sum(maxIvalues >= B(i,1))/nsims);
    end
    
    % Plotting the proportions of simulations in which thresholds for the maximum number of people infected simultaneously are reached
    plot(B(:,1), B(:,2), 'Color', colour(k));
    % Plotting the expected probabilities of major outbreak
    yline(max(0, 1-(1/Rvalues(k))), 'Color', colour(k),'LineStyle','--');

    hold on

end

hold off

% Adding axis labels and a legend
xlabel('Number of individuals infected simultaneously (M)')
ylabel('Proportion of simulations in which threshold M was reached')
legend({'R_0 = 0.8','','R_0 = 1.01','','R_0 = 1.1','','R_0 = 1.3','','R_0 = 1.6','','R_0 = 2','' }, 'Location','northeast')


figure(2)
% Creating a histogram of the maximum number of people infected simultaneously, normalised to give the proportions of simulations where the
% maximum number of people infected is in each bin
histogram(maxIvalues, 25,'Normalization','probability')
xlim([0 250])
xlabel('Number of individuals infected simultaneously')
ylabel('Normalised Frequency Density')
legend('$P(major \; outbreak)=0.5$','Interpreter','latex')





