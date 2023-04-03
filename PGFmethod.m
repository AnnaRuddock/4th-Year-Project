function F = PGFmethod(q)

% Inputting the NGM (from the data from an influenza outbreak in Mexico)
K = [1.41, 0.34; 0.35, 0.87];

% Writing the 2 equations that q(1) and q(2) must satisfy
F(1) = q(1) - 1/(K(1,1)*(1-q(1)) + K(1,2)*(1-q(2))+1);
F(2) = q(2) - 1/(K(2,1)*(1-q(1)) + K(2,2)*(1-q(2))+1);

end



