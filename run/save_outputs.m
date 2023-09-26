function save_outputs(dataset, parcellation)
%% save all the important outputs for easily loading in plotting scripts

%dataset: 'camcan', 'hcp', 'ukb'
%parcellation: '300', 'gm_voxel'

tic

[title_text, ~, data_dir, all_ages, all_subs, correct_T, Neff, thresh, norm_thresh] = feval(dataset, parcellation); 
[psegs, Lams, ages, subs]=save_outputs_analysis(data_dir,all_subs,correct_T,all_ages,Neff,norm_thresh);

toc
[pred_psegs]=numerics_pseg(Neff,thresh,Lams);


save(strcat('output_',parcellation ,'.mat'),'subs','ages','psegs','Lams','pred_psegs')
