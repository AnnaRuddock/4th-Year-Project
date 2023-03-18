P = [6.4;0.019;4;0.16*2.2;0.72;0.25;6.1;0.3;0.4744];

% Increasing axes font size

set(0,'defaultaxesfontsize',16)

% Looping over all parameters in the model
for j = 1:9
    % Defining a step size proportional to the size of parameter being
    % considered
    inc = 10^(floor(log10(P(j)))-2);
    % Creating a vector to store probabilities of major outbreak
    R = zeros(round(2*P(j)/inc), 2);
    % Creating a vector of times at which to undertake calculations
    G = inc:inc:(inc*round(2*P(j)/inc));

    figure (j)
    % Plotting the expected value of the parameter
    xline(P(j), '--', 'Color', [0.5 0 0.5], 'LineWidth', 2)
    hold on
    for i = 1:round(2*P(j)/inc)
        Pmod = P;
        % Modifying the parameter
        Pmod(j) = inc*i;
        % Using the PGF method code to find the probability of no major
        % outbreak
        fun = @(r)PGFmethodebola(r, Pmod);
        x0 = [0,0];
        % Altering the fsolve options to suppress the display
        options = optimoptions('fsolve','Display','none');
        x = fsolve(fun,x0, options);
        % Subtracting from 1 to find the probability of major outbreak
        x = 1-x;
        R(i,:) = x;
    end
    figure(j)
    % Plotting probability of major outbreak starting with 1 infected
    % person in hospital
    plot(G, R(:,1))
    % Plotting probability of major outbreak starting with 1 infected
    % person in the community
    plot(G, R(:,2))

    % Changing by 10% and plotting the new probabilities of major outbreak
    Pinc = P;
    Pdec = P;
    % increasing and decreasing the parameter value by 10%
    Pinc(j) = P(j)*1.1;
    Pdec(j) = P(j)*0.9;

    % Finding the probabilities of no major outbreak in both cases
    fun = @(r)PGFmethodebola(r, Pinc);
    x0 = [0,0];
    options = optimoptions('fsolve','Display','none');
    x = fsolve(fun,x0, options);

    fun2 = @(r)PGFmethodebola(r, Pdec);
    y0 = [0,0];
    options = optimoptions('fsolve','Display','none');
    y = fsolve(fun2,y0, options);

    if x(1)>y(1)
        xline(Pinc(j), '--', 'Color', [34 139 87]/255)
    else
        xline(Pdec(j), '--', 'Color', [34 139 87]/255)
    end
    %if x(2)>y(2)
    %    scatter(Pinc(j), 1-x(2), 25, 'k','filled')
    %else
    %    scatter(Pdec(j), 1-y(2), 25, 'k','filled')
    %end

    hold off
    % Adding y-axis label and legend
    ylabel('Probability')
    %legend({' Expected Parameter Value', 'Probability of major outbreak starting with 1 infected person in hospital', 'Probability of major outbreak starting with 1 infected person in the community'}, Location="northoutside")
    ylim([0 0.8])
   
end

% Adding x-axis labels

figure(1)
xlabel('Number of household contacts, N')
    
figure(2)
xlabel('Hospital transmission rate multiplier, \alpha')

figure(3)
xlabel('Hospital contact multiplier, \beta')

figure(4)
xlabel('Transmission probability, q')

figure(5)
xlabel('Proportion of infected people hospitalised, h')
xlim([0 1])

figure(6)
xlabel('Average number of visitors infected by \newline a hospitalised person, \lambda_h')

figure(7)
xlabel('Average number of people infected at an \newline insecure burial, \phi')

figure(8)
xlabel('Survival probability, \mu')

figure(9)
xlabel('Proportion of burials that are secure, \gamma')
xlim([0 1])
    

