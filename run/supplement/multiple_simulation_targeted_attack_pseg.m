function multiple_simulation_targeted_attack_pseg()


%age 51
outdir = 'output/random/6025360/Lam1.7/';
[psegs, ks] = readin_simulation(outdir);
l = scatter(ks, psegs, 'filled');
hold on

outdir = 'output/targeted/decreasing/density/6025360/Lam1.7/';
[psegs, ks] = readin_simulation(outdir);
l2 = scatter(ks, psegs, 'filled');
hold on

outdir = 'output/targeted/decreasing/length/6025360/Lam1.7/'; 
[psegs, ks] = readin_simulation(outdir);
l3 = scatter(ks, psegs, 'filled');
hold on

%age 65
outdir = 'output/targeted/decreasing/glut4/6025360/Lam1.7/'; 
[psegs, ks] = readin_simulation(outdir);
l4 = scatter(ks, psegs, 'filled');
hold on


hAxis=gca;
hAxis.FontSize = 20;

title('Ising model simulations')
ylabel('{\it P}_{seg}')
xlabel("average degree")

leg = legend([l l2 l3 l4],{'random', 'density','length', 'glut4'}, 'Location', 'southwest');
title(leg,'remove')
