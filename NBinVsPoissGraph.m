x = 0:50;
plot(x,nbinpdf(x,1.024,.3125), ...
     x,poisspdf(x,2.2528));
legend({'Negative Binomial' 'Poisson'})
xlabel('Number of secondary cases')
ylabel('Probability')