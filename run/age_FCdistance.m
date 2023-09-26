function age_FCdistance(dataset, parcellation)
%%check out how functional connections are degraded with respect to age

%dataset: 'camcan', 'hcp', 'ukb'
%parcellation: '300', 'gm_voxel'


dists = load('Seitzman_distance.mat').distances;

short_indi = dists>0 & dists<56;	%56mm delineates first quadrant of distances between regions
long_indi = dists>102;			%106mm delineates third quandrant of distances between regions
all_indi = dists>0;


[title_text, Edges, data_dir, all_ages, all_subs, correct_T, ~, ~, ~] = feval(dataset, parcellation); 

TOTAL_SUBS = size(all_ages,1);

short_corrs = NaN(1, TOTAL_SUBS);
long_corrs = NaN(1, TOTAL_SUBS);
all_corrs = NaN(1, TOTAL_SUBS);
sub_ages = NaN(1, TOTAL_SUBS);


for s=1:TOTAL_SUBS
	sub = all_subs{s};

	filename = get_filename(data_dir, sub);
    
	if isfile(filename)
		if contains(filename, '.mat')
			raw_data = load(filename).raw_data;
		else
			raw_data = readmatrix(filename);
		end	
		
		sz=size(raw_data);
		T=sz(1);	%check to make sure right size!!
		num_regions=sz(2);

		if T == correct_T 
		    R = corrcoef(raw_data);
		    
		    adj_R = R(short_indi);
		    short_cor = mean(adj_R(adj_R>0));
		    short_corrs(s) = short_cor;

		    adj_R = R(long_indi);
		    long_cor = mean(adj_R(adj_R>0));
		    long_corrs(s) = long_cor;

		    adj_R = R(all_indi);
		    all_cor = mean(adj_R(adj_R>0));
		    all_corrs(s) = all_cor;

		    sub_ages(s) = all_ages(s);
		end

	end
end


inds = ~isnan(short_corrs);
sub_ages=sub_ages(inds);
[n_exp, data, rho, pval] = plot_pseg(Edges,short_corrs(inds),sub_ages','b', 'shortest');
[n_exp, data, rho, pval] = plot_pseg(Edges,long_corrs(inds),sub_ages','r', 'longest');
[n_exp, data, rho, pval] = plot_pseg(Edges,all_corrs(inds),sub_ages','k', 'all');



hAxis=gca;
title(title_text)
ylim([0.18 0.35])
ylabel('average correlation')
xlabel('age in years')
xlim([min(Edges)-1 max(Edges)+1])
legend('location', 'best', 'FontSize', 14);

