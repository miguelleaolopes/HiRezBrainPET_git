import subprocess
import os
import re
import csv

def run_script_and_capture_output(image_path, config_path, additional_args):
    command = f"python computeDerenzoValleyToPeak.py -f {image_path} -c {config_path} " + additional_args
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return result.stdout

def extract_vpr_and_resolvability(output):
    lines = output.split('\n')
    avg_vpr_line = next((line for line in lines if line.startswith('Average VPR')), None)
    resolvability_line = next((line for line in lines if line.startswith('Resolvability [%]')), None)
    return avg_vpr_line, resolvability_line

def remove_ansi_codes(text):
    ansi_escape = re.compile(r'\x1B[@-_][0-?]*[ -/]*[@-~]')
    return ansi_escape.sub('', text)

def extract_iteration_from_filename(filename):
    # Assuming filename format is always like '<name>_it<i>.hdr'
    parts = filename.split('_it')
    if len(parts) > 1:
        iteration_part = parts[1]
        iteration_number = iteration_part.split('.')[0]  # Get the part before '.hdr'
        return iteration_number
    return 'Unknown'

def process_images(folder_path, config_path, output_file, metric='parabola', z_range_min='20', z_range_max='30'):
    additional_args = f"-m {metric} -z {z_range_min} {z_range_max}"
    headers = ['Iteration', 'VPR 1.1', 'VPR 1.0', 'VPR 0.95', 'VPR 0.9', 'VPR 0.85', 'VPR 0.8', 'Resol 1.1', 'Resol 1.0', 'Resol 0.95', 'Resol 0.9', 'Resol 0.85', 'Resol 0.8']
    results = []

    if output_file.endswith('.txt'):
        with open(output_file, 'w') as file:
            for image in os.listdir(folder_path):
                if image.endswith('.hdr') and not image.endswith('sensitivity.hdr'):
                    image_path = os.path.join(folder_path, image)
                    output = run_script_and_capture_output(image_path, config_path, additional_args)
                    avg_vpr, resolvability = extract_vpr_and_resolvability(output)
                    if avg_vpr and resolvability:
                        avg_vpr = remove_ansi_codes(avg_vpr)
                        resolvability = remove_ansi_codes(resolvability)
                        file.write(f'{image}\n{avg_vpr}\n{resolvability}\n\n')
    if output_file.endswith('.csv'):
        for image in os.listdir(folder_path):
            if image.endswith('.hdr') and not image.endswith('sensitivity.hdr'):
                iteration = extract_iteration_from_filename(image)
                image_path = os.path.join(folder_path, image)
                output = run_script_and_capture_output(image_path, config_path, additional_args)
                avg_vpr, resolvability = extract_vpr_and_resolvability(output)

                if avg_vpr and resolvability:
                    avg_vpr = remove_ansi_codes(avg_vpr).split(':')[1].split()
                    resolvability = remove_ansi_codes(resolvability).split(':')[1].split()
                    row = [iteration] + avg_vpr + resolvability
                    results.append(row)

        with open(output_file, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(headers)
            for row in results:
                writer.writerow(row)
    
    return output_file

# Example usage
# output_5test = process_images('example/5test_rs', 'example/config.json', 'example/output_5test_rs.csv')
# output_4test = process_images('example/4test_rs', 'example/config.json', 'example/output_4test_rs.csv')

# path_1 = '../../Results/rsmall_z2s/vs015_MLEM_projs/'
# path_2 = '../HiRezBrainPET_git/BPET_castor_results/Results/rsmall_z2s/vs015_opts_dri/D95/'
# path_3 = '../HiRezBrainPET_git/BPET_castor_results/Results/rsmall_z2s/vs015_opts_dri/'
# folder_1 = ['jos','dri','inc','sid']
# folder_2 = ['b01', 'b03', 'b05', 'b08', 'b1', 'b10', 'b100']
# folder_3 = ['MLEM', 'OSL_F/b1', 'OSL_P/b1', 'D95/b1']
# name_1 = 'vs015_MLEM_projs'
# name_2 = 'vs015_opts_dri_D95_betas'
# name_3 = 'vs015_opts_dri'
# outpath = './results'

# path = path_1
# folder = folder_1
# name = name_1
# for f in folder:
#     output_path = outpath+'/output_'+name+'_parabola_z2030.csv'
#     output_file = process_images(path+f, 'example/config.json', output_path)
#     print("Output path: ",output_path)
#     print("")

# path_4 = '../HiRezBrainPET_git/BPET_castor_results/Results/rsmall_z2s/vs015_MLEM_jos_ittest/'
# folder_4 = ['1test', '2test', '3test', '4test', '5test', '6test', '7test', '8test']
# name_4 = 'vs015_MLEM_jos_ittest'
# subsets = [50,25,20,10,5,2,1,[50,10]]


# for f in folder_4:
#     output_path = path_4+f+'/output_'+name_4+'.csv'
#     output_file = process_images(path_4+f, 'example/config.json', output_path)
#     print("Output path: ",output_path)
#     print("")