function multiple_age_pseg(which)
%%plot all 3 datasets' aging trends

%which: 'lambda', 'pseg'

parcellation = '300';

datasets = {'camcan', 'ukb', 'hcp'};
pretty_datasets = {'CamCAN', 'UKB','HCP'};
colors={'c', 'm', 'g'};

for num = 1:length(datasets)
	dataset = datasets{num};
	[title_text, Edges, data_dir, ~, ~, ~, ~,~,~] = feval(dataset, parcellation); 
	filename = strcat('output_',parcellation ,'.mat');

	load(strcat(dataset,'/',filename),'subs','ages','psegs','Lams','pred_psegs')

	inds=Lams<-1;	%remove unphysical fits where lambda is negative (Lamba = -1 means lambda = 0)
		%coupling strength lambda assumed to always be positive
	disp(['number of unphysical lambda fits: ', num2str(sum(inds))])

	ages(inds)=NaN;%NaNs get removed in plot_pseg
	psegs(inds)=NaN;
	Lams(inds)=NaN;

	if contains(lower(which), 'lam')
		ys=Lams;
		y_text = '\Lambda';
		ylim([-0.5 0.1])
	else
		ys=psegs;
		y_text = '{\it P}_{seg}';
		ylim([0.6 1])
	end

	inds = ~isnan(ys);
	ys=ys(inds);
	ages=ages(inds);

	[n_exp, data, rho, pval] = plot_pseg(Edges,ys,ages, colors{num}, pretty_datasets{num});
	hold on
end

ylabel(y_text);
xlabel('age in years')
xlim([19 91])
legend('location', 'best', 'FontSize', 14);

end

