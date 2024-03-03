function fit_Neff(dataset, parcellation)
%%alternative fitting procedure from identify_Neff() where Neff is fitted per individual rather than for all individuals in the dataset

%dataset: 'camcan', 'hcp', 'ukb'
%parcellation: '300', 'gm_voxel'

tic
[title_text, Edges, data_dir, all_ages, subs, correct_T, ~,~,~] = feval(dataset, parcellation);	

[sub_ages, Neffs, m2s, lams] = fit_Neff_analysis(data_dir,subs,correct_T, all_ages);
toc

inds = lams<0;	%get rid of unphysical fits
disp(['number of unphysical lambda fits: ', num2str(sum(inds))])

lams(inds)=NaN;
Neffs(inds)=NaN;
m2s(inds)=NaN;

xlabel_text = 'effective number of regions{\it N}_{eff}';
ylabel_text = 'connection strength \lambda';


%save(strcat(dataset,'_Neffs_',parcellation ,'.mat'),'subs','Neffs')

plot_Neff(Neffs, lams, xlabel_text, ylabel_text, title_text, nnz(~isnan(Neffs)))

end
