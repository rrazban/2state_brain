function [title_text, age_bins, data_dir, ages, subs, correct_T, Neff_opt, threshold, exact_norm_thresh]=camcan(parcellation)

title_text ='Cambridge Centre for Ageing';

%age_bins=20:5:85;	%youngest age is 18 but very few ppl in 18-20 bin
age_bins=17.5:5:87.5;%there is only one person who is 88 who gets excluded by this binning	
%age_bins=18:1:87;

data_dir = strcat('/shared/datasets/public/camcan/derivatives/parcelled/', parcellation);	%gm_voxel, 300
%note this path must be adjusted to match where the processed data is located on your computer
%processed data is not provided as part of this GitHub directory. Files must be directly accessed from Cambridge Center for Ageing and preprocessed according to the atlas of your choice


info_file = 'CAMCAN_participant_data.tsv';
info = tdfread(info_file);
ages = info.age;
subs = cellstr(info.Observations);

correct_T =241;	%its wmcsf but with first 20 timepoints removed to remove coil warmup fluctuations

if contains(data_dir, 'nosmooth')	%gm_voxel with no 5mm FWHM smoothing
	Neff_opt = 280; 
	threshold = 58;
	exact_norm_thresh=0.209036145;
elseif contains(data_dir, 'voxel')
	Neff_opt = 65;	%voxel	%obtained from identify_Neff.m
	threshold = 19;	%Ising model Pseg	
	exact_norm_thresh = 0.297088569;	%exact:19.31075697/65, data Pseg	
elseif contains(data_dir, '300')
	Neff_opt = 40;
	threshold = 12;
	exact_norm_thresh=0.333868715;
end
