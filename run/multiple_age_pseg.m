function multiple_age_pseg(which)
%%plot all 3 datasets' aging trends

%which: 'lambda', 'pseg'

parcellation = '300';	%can also be 'gm_voxel' but data not available for UKB

datasets = {'camcan', 'ukb', 'hcp'};
pretty_datasets = {'CamCAN', 'UKB','HCP'};
colors={'#4DBEEE', '#7E2F8E', '#77AC30'};

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
		loc = 'northeast';
	elseif contains(lower(which), 'pseg')
		ys=psegs;
		y_text = '{\it P}_{seg}';
		ylim([0.6 1])
		loc = 'northeast';
	elseif contains(lower(which), 'pint')
		ys=1-psegs;
		y_text = '{\it P}_{int}';
		ylim([0 0.4])
		loc = 'southeast';

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
legend('location', loc, 'FontSize', 14);

end

