function plot_synch_distribution(dataset, parcellation)
%%plot an individual's synchrony distribution to check if Ising model with Neff captures the data

%dataset: 'camcan', 'hcp', 'ukb'
%parcellation: '300', 'gm_voxel'


[title_text, ~, data_dir, ages, subs, correct_T, Neff,~,~] = feval(dataset, parcellation);

plot_synch_distribution_analysis(data_dir, subs, correct_T, ages, Neff)

