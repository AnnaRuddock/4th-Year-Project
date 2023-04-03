% Inputting a vector of parameter values
P = [6.4;0.019;4;0.16*2.2;0.72;0.25;6.1;0.3;0.4744];

Ebolavary(P)

function[ebolamodify] = Ebolavary(P)

% Creating a vector of parameter names
Parameter_name = ["Number of household contacts";
    "Hospital transmission rate multiplier";
    "Hospital contact multiplier";
    "Transmission probability";
    "Proportion of infected people hospitalised";
    "Average number of visitors infected by a hospitalised person";
    "Average number of people infected at an insecure burial";
    "Survival probability";
    "Proportion of burials that are secure"];
Parameter_symbol = ["N"; "alpha"; "beta"; "q"; "h"; "lambda_h"; "phi"; "mu"; "gamma"];

% Inputting the original parameter values
Original_parameter_value = P;

% Creating vectors to record outputs
Pmajoroutbreak_h_perc_change = zeros(9,1);
Pmajoroutbreak_c_perc_change = zeros(9,1);
Firsthits0_h = zeros(9,1);
Firsthits0_c = zeros(9,1);

    for j = 1:9
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

        % Same calculation but with the estimated parameters
        fun3 = @(r)PGFmethodebola(r, P);
        z0 = [0,0];
        options = optimoptions('fsolve','Display','none');
        z = fsolve(fun3,z0, options);

        % Finding the minimum probability of major outbreak starting in
        % each treatment class
        minPhospital = max([0, 1-max([x(1),y(1)])]);
        minPcommunity = max([0,1-max([x(2),y(2)])]);

        % Saving these values so they can be used later
        holdx=x;
        holdy=y;

        % Finding the percentage changes in the probability of major
        % outbreak from the expected values
        Pmajoroutbreak_h_perc_change(j) = (minPhospital - (1-z(1)))/(1-z(1))*100;
        Pmajoroutbreak_c_perc_change(j) = (minPcommunity - (1-z(2)))/(1-z(2))*100;

        % Defining step size, and vectors to store probabilities and
        % timesteps
        inc = 10^(floor(log10(P(j)))-2);
        R = zeros(round(2*P(j)/inc), 2);
        G = inc:inc:(inc*round(2*P(j)/inc));
          
        % Finding probability of major outbreak over a range of parameter
        % values, from 0 to 2*Expected value (for more detail read
        % Ebolaplots file)
        for i = 1:round(2*P(j)/inc)
            Pmod = P;
            Pmod(j) = inc*i;
            fun = @(r)PGFmethodebola(r, Pmod);
            x0 = [0,0];
            options = optimoptions('fsolve','Display','none');
            x = fsolve(fun,x0, options);
            x = 1-x;
            R(i,:) = x;
        end

        % Finding the first time at which the probability of major outbreak
        % hits 0, if this occurs

        % Using if statements to find whether increasing or decreasing the
        % parameter leads to a lower probability of major outbreak
        % 
        if min(R(:,1))>0.0001
            Firsthits0_h(j) = "N/A";
        else if holdx(1)>holdy(1)
            % Using the find function to find the first time at which the
            % probability of major outbreak is very close to 0 (it will
            % never actually hit zero)
            Firsthits0_h(j) = find(R(:,1)<0.0001, 1)*inc;
        elseif holdx(1)<holdy(1)
            Firsthits0_h(j) = find(R(:,1)<0.0001, 1, 'last')*inc;
        end

        % Same calculations but for starting in the community rather than
        % in hospital
        if min(R(:,2))>0.0001
            Firsthits0_c(j) = "N/A"
        elseif holdx(2) > holdy(2)
            Firsthits0_c(j) = find(R(:,2)<0.0001, 1)*inc;
        elseif holdx(2) < holdy(2)
            Firsthits0_c(j) = find(R(:,2)<0.0001, 1, 'last')*inc;
        end

        end
          
    end

% Collating all data into a table
ebolamodify = table(Parameter_name, Parameter_symbol, Original_parameter_value, Pmajoroutbreak_h_perc_change, Pmajoroutbreak_c_perc_change, Firsthits0_h, Firsthits0_c);

% Sorting by the percentage change in the probability of major outbreak
% starting with 1 infected person being treated in the community
ebolamodify = sortrows(ebolamodify, 5)

end
