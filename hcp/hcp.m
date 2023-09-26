function [title_text, age_bins, data_dir, ages, subs, correct_T, Neff_opt, threshold, exact_norm_thresh]=hcp(parcellation)

title_text ='Human Connectome Project';
age_bins = 35:5:90;
%age_bins = 35:2:90;

data_dir = strcat('/shared/datasets/public/hcp-a/derivatives/parcelled/', parcellation);	%gm_voxel, Rostam_300/
%300/ by Botond requires extra readin preprocessing	%ask to rename??

info_file = 'info_hcp.txt';
info = tdfread(info_file);
ages = info.age;
subs = cellstr(info.src_subject_id);

correct_T = 1912;

if contains(data_dir, 'nosmooth')	%gm_voxel with no 5mm FWHM smoothing
	Neff_opt = 950; 
	threshold = 146;
	exact_norm_thresh = 0.154258373;
elseif contains(data_dir, 'voxel')
	Neff_opt = 125; 
	threshold = 31;	%Ising model Pseg	
	exact_norm_thresh = 0.254368421;	%exact:31.7960526/125, data Pseg	
elseif contains(data_dir, '300')
	Neff_opt = 40;
	threshold = 12;
	exact_norm_thresh=0.333868715;
end
