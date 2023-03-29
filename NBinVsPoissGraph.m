set(0,'defaultaxesfontsize',12)

x = 0:15;
plot(x,nbinpdf(x,1.024,.3125), ...
     x,poisspdf(x,2.2528), x,geopdf(x,0.3074));
legend({'Negative Binomial (1.024, 0.3125) Distribution' 'Poisson (2.2528) Distribution', 'Geometric(0.3074) Distribution'})
xlabel('Number of secondary cases')
ylabel('Probability')
