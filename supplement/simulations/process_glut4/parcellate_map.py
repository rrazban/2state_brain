"""

Code to process gene expression data into 
atlas of choice. Gene expression data is 
manually downloaded from AHBA
https://www.meduniwien.ac.at/neuroimaging/mRNA.html

"""

import sys
import numpy as np
import nibabel as nib
import nilearn


def get_atlas_voxel(filename):
	roi_img = nib.load(filename)
	roi_data = roi_img.get_data()
	return roi_data

if __name__ == '__main__':
	cort_roi_file = 'atlas/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz'
	subcort_roi_file = 'atlas/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz'
	map_file = '6517/6517_mirr_mRNA.nii'

	cort_Har_roi = get_atlas_voxel(cort_roi_file)
	subcort_Har_roi = get_atlas_voxel(subcort_roi_file)
	gene_map = get_atlas_voxel(map_file)

	max_x, max_y, max_z = gene_map.shape

	rois = np.zeros(71)
	cort_offset=22


	special_subcortical_ignore = [0, 1, 2, 8, 12, 13]	#either white matter or includes cortex (broadly defined) #also brainstem here cuz no white matter
	white_matter = [1,12]
	cortex = [2,13]
	ignore = [0,1,8,12]	#0 is background, 8 is brainstem	
	for x in range(max_x):
		map_x = max_x-1 - x	#for some reason has diff convention only for x	
		for y in range(max_y):
			for z in range(max_z):
				gene_density = gene_map[map_x, y, z]
				roi = int(subcort_Har_roi[x, y, z])	#for some reason a float	#cort_Har_roi does not have this problem
				if roi in cortex:
					cortex_roi = cort_Har_roi[x,y,z]+cort_offset
					if cortex_roi!=cort_offset:#cooresponds to background		#i should check inverse too that subcort includes everything
						rois[cortex_roi]+=gene_density
				elif roi in ignore:
					continue
				else:
					rois[roi]+=gene_density


	d_out = {}
	num=0
	for i, val in enumerate(rois):
		if val!=0:
			d_out[num] = val
			num+=1

#	d_out = {i:rois[i] for i in range(len(rois))}
	sorted_x = sorted(d_out.items(), key=lambda kv: kv[1], reverse=True)
	print(sorted_x)




