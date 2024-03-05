function simulation_age_pseg()


outdir = 'output/random/6025360_20250_2_0_density/Lam1.7/';
[psegs, ks] = readin_simulation(outdir);

Edges=min(ks):1:max(ks)+1;
hl = scatter(ks, psegs,50,[0.8500 0.3250 0.0980], 'filled');
hold on
[n_exp, l_exp] = plot_pseg(Edges,psegs,ks,'m','');


title('Ising model simulations')
ylabel('{\it P}_{seg}')
xlabel('average degree')
[rho,pval]=corr(ks,psegs','Type','Spearman')
legend([hl l_exp],{strcat('\rho=', num2str(rho,'%.3f'),' (', num2str(pval,2), ') \newlineN=', string(n_exp)), 'binned data'}, 'Location', 'southwest')
