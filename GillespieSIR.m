clc
clear
close all

Rvalues = [0.8, 1.01, 1.1, 1.3, 1.6, 2];
colour = ["r" "g" "b" "c" "m" "k"];

for k = 1:6

    u = 1;
    N = 1000;
    beta = (Rvalues(k)*u)/N;
    nsims = 10000;
    maxIvalues = zeros(nsims,1);
    
    for j = 1:nsims

        if mod(j,100000) == 0
            j
        end

        S = 999;
        I = 1;
        maxI = I;
        R = 0;
        while I>0
            lambda = ((beta*S*I))+u*I;

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
        
        maxIvalues(j) = maxI;
    
    end
    
    B = zeros(251,2);
    B(:, 1) = 0:250;
    
    for i = 1:251
        B(i,2) = (sum(maxIvalues >= B(i,1))/nsims);
    end
    
    plot(B(:,1), B(:,2), 'Color', colour(k));
    yline(max(0, 1-(1/Rvalues(k))), 'Color', colour(k),'LineStyle','--');

    hold on

end

xlabel('Number of individuals infected simultaneously (M)')
ylabel('Proportion of simulations in which threshold M was reached')
legend({'R_0 = 0.8','','R_0 = 1.01','','R_0 = 1.1','','R_0 = 1.3','','R_0 = 1.6','','R_0 = 2','' }, 'Location','northeast')









