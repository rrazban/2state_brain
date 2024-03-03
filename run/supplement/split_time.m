function split_time()
%% perform Ising model analysis for 5 chunks of time within fMRI scan 
%investigate temporal robustness of pseg


parcellation='300';

datasets = {'camcan', 'ukb', 'hcp'};
pretty_datasets = {'CamCAN', 'UKB','HCP'};
colors={'#4DBEEE', '#7E2F8E', '#77AC30'};


tic
for num = 1:length(datasets)
	dataset = datasets{num};

	if contains(dataset, 'hcp')
		parcellation='Rostam_300'
	end
	[title_text, Edges, data_dir, all_ages, all_subs, correct_T, Neff, thresh, norm_thresh] = feval(dataset, parcellation); 
	[psegs, Lams, ages, subs]=split_time_analysis(data_dir,all_subs,correct_T,all_ages,Neff,norm_thresh);



	inds=Lams<-1;	%remove unphysical fits where lambda is negative (Lamba = -1 means lambda = 0)
%	ages(inds)=NaN;	%cant do this cuz wrong size
	psegs(inds)=NaN;
	Lams(inds)=NaN;

	ys = std(psegs, 1);

	[n_exp, data, rho, pval] = plot_pseg(Edges,ys,ages, colors{num}, pretty_datasets{num});

end
toc
y_text = 'standard deviation of {\it P}_{seg}';
ylabel(y_text);
xlabel('age in years')
xlim([19 91])
ylim([0 0.105])
legend('location', 'southwest', 'FontSize', 14);

end

