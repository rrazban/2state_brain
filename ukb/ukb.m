function [title_text, age_bins, data_dir, ages, subs, correct_T, Neff_opt, threshold, exact_norm_thresh]=ukb(parcellation)

title_text ='UK Biobank';
age_bins = 45:5:80;	

data_dir = strcat('/shared/datasets/public/ukb/derivatives/parcelled/', parcellation,'/001hz/');	%300

info_file = 'phenotypes.csv';
info = tdfread(info_file,',');
ages = info.age;
subs = cellstr(num2str(info.eid));

correct_T = 490;


if contains(data_dir, '300')
	Neff_opt = 30;
	threshold = 10;
	exact_norm_thresh=0.357089829;
end
