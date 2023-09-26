function filename = get_filename(data_dir, sub)
	 %can accomodate bolus filename as well
	if contains(data_dir, 'voxel')
		filename = strcat(data_dir,'/sub-',sub,'.mat');	%voxel data in .mat cuz too big to load .csv
	else
		filename = strcat(data_dir,'/sub-',sub,'.csv'); 
	%	filename = strcat(data_dir,'/',sub,'.csv');	%for UKB 
	end
end
