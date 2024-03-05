function multiple_simulation_age_pseg()


%age 51
outdir = 'output/random/6025360/Lam1.7/';
[psegs, ks] = readin_simulation(outdir);
l = scatter(ks, psegs, 'filled');
hold on

%age 57
outdir = 'output/random/4712851/Lam1.7/';
[psegs, ks] = readin_simulation(outdir);
l2 = scatter(ks, psegs, 'filled');
hold on

%age 61
outdir = 'output/random/3081886/Lam1.7/';
[psegs, ks] = readin_simulation(outdir);
l3 = scatter(ks, psegs, 'filled');
hold on

%age 65
outdir = 'output/random/1471888/Lam1.7/';
[psegs, ks] = readin_simulation(outdir);
l4 = scatter(ks, psegs, 'filled');
hold on

%age 72
outdir = 'output/random/4380337/Lam1.7/';
[psegs, ks] = readin_simulation(outdir);
l5 = scatter(ks, psegs, 'filled');
hold on

%age 74
outdir = 'output/random/1003054/Lam1.7/';
[psegs, ks] = readin_simulation(outdir);
hl = scatter(ks, psegs, 'filled');
hold on


hAxis=gca;
hAxis.FontSize = 20;

title('Ising model simulations')
ylabel('{\it P}_{seg}')
xlabel("average degree")

leg = legend([l l2 l3 l4 l5 hl ],{'51', '57', '61', '65', '72', '74'}, 'Location', 'southwest');
title(leg,'age')
