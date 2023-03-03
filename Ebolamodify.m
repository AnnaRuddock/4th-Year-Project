function[ebolamodify] = Ebolamodify(P)

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

Original_parameter_value = P;

%}
Pmajoroutbreak_h_perc_change = zeros(9,1);
Pmajoroutbreak_c_perc_change = zeros(9,1);
Firsthits0_h = zeros(9,1);
Firsthits0_c = zeros(9,1);

    for j = 1:9
        Pinc = P;
        Pdec = P;
        Pinc(j) = P(j)*1.1;
        Pdec(j) = P(j)*0.9;

        fun = @(r)PGFmethodebola(r, Pinc);
        x0 = [0,0];
        options = optimoptions('fsolve','Display','none');
        x = fsolve(fun,x0, options);

        fun2 = @(r)PGFmethodebola(r, Pdec);
        y0 = [0,0];
        options = optimoptions('fsolve','Display','none');
        y = fsolve(fun2,y0, options);

        fun3 = @(r)PGFmethodebola(r, P);
        z0 = [0,0];
        options = optimoptions('fsolve','Display','none');
        z = fsolve(fun3,z0, options);

        minPhospital = max([0, 1-max([x(1),y(1)])]);
        minPcommunity = max([0,1-max([x(2),y(2)])]);

        holdx=x;
        holdy=y;

        Pmajoroutbreak_h_perc_change(j) = (minPhospital - (1-z(1)))/(1-z(1))*100;
        Pmajoroutbreak_c_perc_change(j) = (minPcommunity - (1-z(2)))/(1-z(2))*100;

        inc = 10^(floor(log10(P(j)))-2);
        R = zeros(round(2*P(j)/inc), 2);
        G = inc:inc:(inc*round(2*P(j)/inc));
          
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
        R;
        %min(R(:,1))
        %find(R(:,1)<0.001)



        if min(R(:,1))>0.0001
            Firsthits0_h(j) = "N/A";
        else if holdx(1)>holdy(1)
            Firsthits0_h(j) = find(R(:,1)<0.0001, 1)*inc;
        elseif holdx(1)<holdy(1)
            Firsthits0_h(j) = find(R(:,1)<0.0001, 1, 'last')*inc;
        end

        if min(R(:,2))>0.0001
            Firsthits0_c(j) = "N/A"
        elseif holdx(2) > holdy(2)
            Firsthits0_c(j) = find(R(:,2)<0.0001, 1)*inc;
        elseif holdx(2) < holdy(2)
            Firsthits0_c(j) = find(R(:,2)<0.0001, 1, 'last')*inc;
        end
          
    end
    
ebolamodify = table(Parameter_name, Parameter_symbol, Original_parameter_value, Pmajoroutbreak_h_perc_change, Pmajoroutbreak_c_perc_change, Firsthits0_h, Firsthits0_c);

ebolamodify = sortrows(ebolamodify, 5)

end