function [psegs,avg_degrees]=readin_simulations(out_dir)

base = '/shared/home/rostam/int_seg/supplement/simulations';
data_dir = strcat(base,'/', out_dir);


info_file = strcat(data_dir, '/info.csv');
info = tdfread(info_file,',');
ks = info.avg_degree;

subs=compose('%d', info.iteration);    %works when have single digit and two digits

correct_T = 2501; 
Neff = 64;  %simply equal to number of nodes with edges in parcellation 
norm_thresh = 0.299;


[psegs, Lams, avg_degrees, sub_ids] = simulation_analysis(data_dir,subs,correct_T,ks,Neff,norm_thresh);

