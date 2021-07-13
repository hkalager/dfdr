function pi_0hat=compute_pi_0hat(pvalues, lambda)

W=sum(pvalues>lambda);

nb_strats=size(pvalues,1);

pi_0hat=min(1, W/((1-lambda)*nb_strats));
