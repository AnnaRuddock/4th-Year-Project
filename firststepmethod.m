function F = firststepmethod(q)

% inputting the NGM matrix (taken from a model of influenza in Mexico)
K = [1.41, 0.34; 0.35, 0.87];

% equations that q(1), q(2) must satisfy
F(1) = ((K(1,1))*(q(1))^2 + K(1,2)*(q(1))*(q(2)) + 1)/(K(1,1)+K(1,2)+1) - q(1);
F(2) = ((K(2,1))*(q(1))*(q(2))+ (K(2,2)*(q(2))^2 + 1))/(K(2,1) + K(2,2) + 1) - q(2);

end
