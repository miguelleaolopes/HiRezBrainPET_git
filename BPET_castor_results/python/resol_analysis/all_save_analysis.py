from save_derenzo_valleytopeak import process_images
from analysis_plot_VPR_Resol import plot_vpr_resol, plot_vpr_resol_it_subsets

path_ittest = '../../Results/rsmall_z2s/vs015_MLEM_jos_ittest/'
folder_ittest = ['1test', '2test', '3test', '4test', '5test', '6test', '7test', '8test']
name_ittest = 'vs015_MLEM_jos_ittest'
labels_ittest = ['10:50', '20:25', '25:20', '50:10', '100:5', '250:2', '500:1', '2:50,40:10']
subsets = [50,25,20,10,5,2,1,[50,10]]
subsets = [50,25,20,10,5,2,1,10]

path_proj = '../../Results/rsmall_z2s/vs015_MLEM_projs/'
path_d95f = '../../Results/rsmall_z2s/vs015_opts_dri/D95/'
path_oslp = '../../Results/rsmall_z2s/vs015_opts_dri/OSL_P/'
path_oslf = '../../Results/rsmall_z2s/vs015_opts_dri/OSL_F/'
path_opts = '../../Results/rsmall_z2s/vs015_opts_dri/'
path_conv = '../../Results/rsmall_z2s/vs015_MLEM_dri_g12_3s/'
path_geom = '../../Results/'
# path_4 = '../../Results/rsmall_z2s/vs015_MLEM_dri_g12_3s/'

folder_proj = ['jos','dri','inc','sid']
folder_d95f = ['b01', 'b03', 'b05', 'b08', 'b1', 'b10', 'b100']
folder_oslp = ['b01', 'b05', 'b1', 'b10', 'b100']
folder_oslf = ['5test', '7test', '2test', '1test', 'b1']
folder_opts = ['MLEM', 'OSL_F/b1', 'OSL_P/b1', 'D95/b100']
folder_conv = ['1test', '2test', '3test', '4test', '5test', '6test', '7test', '8test', '9test', '10test', '11test']
folder_geoms = ['rsmall_orig','rsmall_z2s', 'rsmall_big_z4s', 'rsmall_big', 'rmid_orig','rmid_z2s', 'rmid_big_z4s', 'rmid_big', 'rbig_orig','rbig_z2s']
folder_geoms_2 = '/vs015_geoms/1test'
folder_geom = [elem + folder_geoms_2 for elem in folder_geoms]

out_proj = folder_proj
out_d95f = folder_d95f
out_oslp = folder_oslp
out_oslf = ['b01', 'b1', 'b10', 'b100', 'b1-2']
out_opts = ['MLEM', 'OSL_F_b1', 'OSL_P_b1', 'D95_b1']
out_conv = folder_conv
out_geom = folder_geoms
name_proj = 'vs015_MLEM_projs_parabola_z2030'
name_d95f = 'vs015_opts_dri_D95_betas_parabola_z2030'
name_oslp = 'vs015_opts_dri_OSL_P_betas_parabola_z2030'
name_oslf = 'vs015_opts_dri_OSL_F_betas_parabola_z2030'
name_opts = 'vs015_opts_dri_parabola_z2030'
name_conv = 'vs015_convs_MLEM_dri_g12_3s_parabola_z2030'
name_geom = 'vs015_geoms_parabola_z2030'
outpath = './results/'

rod_sizes = [1.1, 1.0, 0.95, 0.9, 0.85, 0.8]  # Rod sizes
rod_sizes_1 = [1.1, 1.0, 0.95]  # Less rod sizes
rod_sizes_2 = [0.9, 0.85, 0.8]  # Less rod sizes
rod_sizes_3 = [0.9]  # Less rod sizes

out = out_geom
path = path_geom
folder = folder_geom
name = name_geom
output_path = []
for i, f in enumerate(folder):
    output_path.append(outpath+"output_"+out[i]+'_'+name+'.csv')


# ----- saving csv files -----
for i, f in enumerate(folder):
    output_path = outpath+"output_"+out[i]+'_'+name+'.csv'
    print("Begin Output path: ",output_path)
    output_file = process_images(path+f, 'example/config.json', output_path, 'parabola', '20', '30')
    print("End Output path: ",output_path)
    print("")

# for f in folder_ittest:
#     output_path = outpath+'output_'+f+'_'+name_ittest+'.csv'
#     output_file = process_images(path_ittest+f, 'example/config.json', output_path)
#     print("Output path: ",output_path)
#     print("")

# ----- plotting -----

output_paths = [[], [], []]
paths = [path_proj, path_d95f, path_oslp, path_opts]
folders = [folder_proj, folder_d95f, path_oslp, folder_opts]
names = [name_proj, name_d95f, path_oslp, name_opts]

# for i, (path, folder, name) in enumerate(zip(paths, folders, names)):
#     for f in folder:
#         output_paths[i].append(outpath+out[i]+name + '.csv')

# -----------------------------
# PROJECTORS
out = out_proj
path = path_proj
folder = folder_proj
name = name_proj
output_path = []
for i, f in enumerate(folder):
    output_path.append(outpath+"output_"+out[i]+'_'+name+'.csv')
# plot_vpr_resol(output_path, folder, rod_sizes)
# plot_vpr_resol(output_path, folder, rod_sizes_1,figsize=(8, 9),fontsize=12,legend_loc=["upper right", "lower right"], save_bool=True, save_path=outpath+name+'_1.png')
# plot_vpr_resol(output_path, folder, rod_sizes_2,figsize=(8, 9),fontsize=12,legend_loc=["lower left", "upper left"], save_bool=True, save_path=outpath+name+'_2.png')
# plot_vpr_resol(output_path, folder, rod_sizes_3)

# -----------------------------
# BETAS - D95F
out = out_d95f
path = path_d95f
folder = folder_d95f
name = name_d95f
output_path = []
for i, f in enumerate(folder):
    output_path.append(outpath+"output_"+out[i]+'_'+name+'.csv')
# plot_vpr_resol(output_path, folder, rod_sizes_1,figsize=(8, 9),fontsize=12,legend_loc=["lower left", "lower center"], save_bool=True, save_path=outpath+name+'_1.png')
# plot_vpr_resol(output_path, folder, rod_sizes_2,figsize=(8, 9),fontsize=12,legend_loc=["lower left", "upper left"], save_bool=True, save_path=outpath+name+'_2.png')

# BETAS - OSLP
out = out_oslp
path = path_oslp
folder = folder_oslp
name = name_oslp
output_path = []
for i, f in enumerate(folder):
    output_path.append(outpath+"output_"+out[i]+'_'+name+'.csv')
# plot_vpr_resol(output_path, folder, rod_sizes_1,figsize=(8, 9),fontsize=12,legend_loc=["lower left", "lower center"], save_bool=True, save_path=outpath+name+'_1.png')
# plot_vpr_resol(output_path, folder, rod_sizes_2,figsize=(8, 9),fontsize=12,legend_loc=["lower left", "upper left"], save_bool=True, save_path=outpath+name+'_2.png')

# -----------------------------
# OPTIMIZERS
out = out_opts
path = path_opts
folder = folder_opts
name = name_opts
output_path = []
for i, f in enumerate(folder):
    output_path.append(outpath+"output_"+out[i]+'_'+name+'.csv')
# plot_vpr_resol(output_path, folder, rod_sizes_1,figsize=(8, 9),fontsize=12,legend_loc=["lower left", "lower center"], save_bool=True, save_path=outpath+name+'_1.png')
# plot_vpr_resol(output_path, folder, rod_sizes_2,figsize=(8, 9),fontsize=12,legend_loc=["lower left", "upper left"], save_bool=True, save_path=outpath+name+'_2.png')


# -----------------------------
# ITERATIONS
output_paths_its = [outpath + 'output_' + f + '_' + name_ittest + '.csv' for f in folder_ittest]
# plot_vpr_resol_it_subsets(output_paths_its, labels_ittest, subsets, rod_sizes,figsize=(8, 9),fontsize=12,legend_loc=["upper right", "upper left"], save_bool=True, save_path=outpath + name_ittest +'.pdf')
# plot_vpr_resol_it_subsets(output_paths_its, labels_ittest, subsets, rod_sizes_1,figsize=(8, 9),fontsize=12,legend_loc=["upper right", "upper left"], save_bool=True, save_path=outpath + name_ittest +'_1.pdf')
# plot_vpr_resol_it_subsets(output_paths_its, labels_ittest, subsets, rod_sizes_2,figsize=(8, 9),fontsize=12,legend_loc=["upper right", "upper left"], save_bool=True, save_path=outpath + name_ittest +'_2.pdf')
# plot_vpr_resol_it_subsets(output_paths_its, labels_ittest, subsets, rod_sizes_3,figsize=(8, 9),fontsize=12,legend_loc=["upper right", "upper left"], save_bool=False, save_path=outpath + name_ittest +'_3.pdf')
