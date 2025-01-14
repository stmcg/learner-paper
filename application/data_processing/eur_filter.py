import os
import pandas as pd
import gzip
import numpy as np


result_table_path = "filtered_table_EUR.txt"
list_df = pd.read_csv("SNP_name_test_only.csv")
snp_names = set(list_df['x'].astype(str).tolist())
# Initialize an empty DataFrame or read the existing one
if os.path.exists(result_table_path):
    result_table = pd.read_csv(result_table_path, sep='\t', index_col=0)
else:
    result_table = pd.DataFrame()
    result_table.to_csv(result_table_path, sep='\t')
    
def process_folder(folder_path, disease_name):
    global result_table
    processed = False
    for root, _, files in os.walk(folder_path):
        for file_name in files:
            if file_name.endswith("auto.txt.gz"):
                file_path = os.path.join(root, file_name)
                if not processed:
                    print(file_path)
                    with gzip.open(file_path, 'rt') as f:
                        # Read in chunks to handle large files
                        df = pd.read_csv(f, sep='\t', usecols=['CHR', 'POS', 'BETA', 'SE'])

                        df['CHR_POS'] = df['CHR'].astype(str) + ':' + df['POS'].astype(int).astype(str)
                        df['Z'] = df['BETA'] / df['SE']
                        temp_dict = {}
                        if disease_name not in result_table.columns:
                            result_table[disease_name] = np.nan
                        for _, row in df.iterrows():
                            chr_pos = f"{row['CHR']}:{int(row['POS'])}"
                            if chr_pos in snp_names:
                                z_score = row['BETA'] / row['SE']
                                if chr_pos not in temp_dict:
                                    temp_dict[chr_pos] = {disease_name: z_score}
                                else:
                                    temp_dict[chr_pos][disease_name] = z_score
                        temp_df = pd.DataFrame.from_dict(temp_dict, orient='index')
                        result_table = result_table.combine_first(temp_df)
                        del df, temp_df    

                    # Set the flag to True after processing the first file
                    processed = True
                if processed:
                    break
        if processed:
            break

def main():
    current_folder = os.path.dirname(os.path.realpath(__file__))
    csv_file_path = os.path.join(current_folder, 'EUR_pheno_name.csv')
    phe_df = pd.read_csv(csv_file_path)
    phe_names = phe_df['x'].tolist()
    target_folder = os.path.join(current_folder, "220_phenotype_download")
    save_counter = 0  # Counter to track how many diseases have been processed
    for dir_name in next(os.walk(target_folder))[1]:
        dir_path = os.path.join(target_folder, dir_name)
        if os.path.isdir(dir_path) and "EUR" in dir_name:
            parts = dir_name.split('.')
            disease_name = parts[-2] if len(parts) > 1 else 'unknown'
            if disease_name in phe_names and disease_name not in result_table.columns:
                process_folder(dir_path, disease_name)
                save_counter += 1
                if save_counter % 10 == 0:  # Save every 10 processed diseases
                    save_path = f"filtered_table_EUR_{save_counter}.txt"  # Adjust the file name here
                    result_table.to_csv(save_path, sep='\t')
                    print(f"Saved result_table after processing {save_counter} diseases to {save_path}.")
    result_table.to_csv(result_table_path, sep='\t')
    print("Final result_table saved.")

if __name__ == "__main__":
    main()