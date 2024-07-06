import pandas as pd

# Load CSV files
path = './../../OneDrive - Universidade de Lisboa/College/MEFT/5Ano/Dissertacao/RPC-BPET/Reconstrucao_PET/CASToR/HiRezBrainPET_git/BPET_castor_results/python/resol_analysis/results/'
df1 = pd.read_csv(path + 'output_1test_vs015_MLEM_jos_ittest.csv')
df2 = pd.read_csv(path + 'output_2test_vs015_MLEM_jos_ittest.csv')
df4 = pd.read_csv(path + 'output_4test_vs015_MLEM_jos_ittest.csv')
df5 = pd.read_csv(path + 'output_5test_vs015_MLEM_jos_ittest.csv')
df6 = pd.read_csv(path + 'output_6test_vs015_MLEM_jos_ittest.csv')
df7 = pd.read_csv(path + 'output_7test_vs015_MLEM_jos_ittest.csv')
df8 = pd.read_csv(path + 'output_8test_vs015_MLEM_jos_ittest.csv')

# Sort DataFrames by 'Iteration'
df1 = df1.sort_values(by='Iteration')
df2 = df2.sort_values(by='Iteration')
df4 = df4.sort_values(by='Iteration')
df5 = df5.sort_values(by='Iteration')
df6 = df6.sort_values(by='Iteration')
df7 = df7.sort_values(by='Iteration')
df8 = df8.sort_values(by='Iteration')

# Function to calculate step size
def calculate_step(total_rows, desired_iterations):
    return total_rows // desired_iterations

# Function to determine the starting index
def starting_index(total_rows, desired_iterations):
    step = calculate_step(total_rows, desired_iterations)
    return step - 1  # Start from the second index

def subsample_replace_iteration(df, desired_iterations):
    subsampled_df = df.iloc[starting_index(len(df), desired_iterations)::calculate_step(len(df), desired_iterations)]
    subsampled_df.loc[:, 'Iteration'] = range(1, 11)
    return subsampled_df

def subtract_dfs(df1, df2):
    diff_df = df1.copy()
    diff_df.set_index('Iteration', inplace=True)
    diff_df.sort_index(inplace=True)
    diff_df = diff_df.subtract(df2)
    return diff_df

# Subsample and replace 'Iteration' column values with 1 to 10
subsampled_df1 = subsample_replace_iteration(df1, 10)
subsampled_df2 = subsample_replace_iteration(df2, 10)
subsampled_df4 = subsample_replace_iteration(df4, 10)
subsampled_df5 = subsample_replace_iteration(df5, 10)
subsampled_df6 = subsample_replace_iteration(df6, 10)
subsampled_df7 = subsample_replace_iteration(df7, 10)
# d8 is different
d8_rows_to_select = [0, 1] + list(range(6, 42, 5))
subsampled_df8 = df8.iloc[d8_rows_to_select]
subsampled_df8.loc[:, 'Iteration'] = range(1, 11)

# Calculate the difference between each DataFrame and subsampled_df1
diff_df1 = subsampled_df1.copy()
diff_df1.set_index('Iteration', inplace=True)
diff_df1.sort_index(inplace=True)

diff_df2 = subtract_dfs(subsampled_df2, diff_df1)
diff_df4 = subtract_dfs(subsampled_df4, diff_df1)
diff_df5 = subtract_dfs(subsampled_df5, diff_df1)
diff_df6 = subtract_dfs(subsampled_df6, diff_df1)
diff_df7 = subtract_dfs(subsampled_df7, diff_df1)
diff_df8 = subtract_dfs(subsampled_df8, diff_df1)

diff_df1 = diff_df1.subtract(diff_df1)

# Write subsampled data and differences to new CSV files
diff_df1.to_csv(path + 'output_1test_vs015_MLEM_jos_ittest_diff.csv')
diff_df2.to_csv(path + 'output_2test_vs015_MLEM_jos_ittest_diff.csv')
diff_df4.to_csv(path + 'output_4test_vs015_MLEM_jos_ittest_diff.csv')
diff_df5.to_csv(path + 'output_5test_vs015_MLEM_jos_ittest_diff.csv')
diff_df6.to_csv(path + 'output_6test_vs015_MLEM_jos_ittest_diff.csv')
diff_df7.to_csv(path + 'output_7test_vs015_MLEM_jos_ittest_diff.csv')
diff_df8.to_csv(path + 'output_8test_vs015_MLEM_jos_ittest_diff.csv')