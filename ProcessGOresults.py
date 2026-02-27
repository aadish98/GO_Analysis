import pandas as pd
import os
import glob
from openpyxl import Workbook

def replace_fbgn_with_symbols(fbgn_ids, fbgn_to_symbol):
    # Handle cases where fbgn_ids is NaN or not a string
    if pd.isna(fbgn_ids) or not isinstance(fbgn_ids, str):
        print("HERERER")
        return 'Unknown'

    # Split the string by commas and convert each id to uppercase
    fbgn_ids_list = fbgn_ids.replace(" ", "").split(',')
    fbgn_ids_list = [x.upper().strip() for x in fbgn_ids_list]

    # Create a dictionary for mapping primary_FBid to current_symbol


    # Map each fbgn_id to the corresponding current_symbol or 'Unknown' if not found
    symbols = [fbgn_to_symbol.get(fbgn_id, 'Unknown') for fbgn_id in fbgn_ids_list]

    # Return the mapped symbols as a comma-separated string
    return ', '.join(symbols)


def process_csv_files(directory_path, lookup_dict, convert_to_name):
    # Process XLSX GO analysis files in the directory (stop writing CSVs)
    fbgn_to_symbol = dict(zip(lookup_dict['primary_FBid'].str.upper(), lookup_dict['current_symbol']))
    xlsx_files = glob.glob(os.path.join(directory_path, '**.xlsx'), recursive=True)

    for xlsx_file in xlsx_files:
        print(f"Processing {xlsx_file}")
        if "GO_Analysis" not in xlsx_file or "summary" in xlsx_file:
            print("Skipped")
            continue  # Skip this file and move to the next one
        # Load the XLSX file
        df = pd.read_excel(xlsx_file)
        if df.shape[0] == 0:
            print("EMpty df")
            continue
        # Check if 'Genes' column exists  
        if 'Genes' in df.columns:
            # Replace FBgn IDs with gene symbols

            df['Genes'] = df['Genes'].apply(
                                        lambda x: replace_fbgn_with_symbols(x, fbgn_to_symbol) if convert_to_name
                                        else x.replace('GN', 'gn')
                                             )

            # Only keep and proceed with significant terms (FDR ≤ 10) if FDR exists
            if 'FDR' in df.columns:
                df = df[df['FDR'] <= 10]
                if df.shape[0] == 0:
                    print("No significant terms (FDR < 10%). Skipping write for", csv_file)
                    continue
                df_sorted = df.sort_values(by='FDR', ascending=True)
            else:
                # If FDR column does not exist, skip writing outputs
                print("No FDR column present. Skipping write for", csv_file)
                continue
            if 'Category' in df_sorted.columns:
            # Group by 'Category' and save each group as a separate sheet in an XLSX file
                output_file = os.path.splitext(xlsx_file)[0] + ".xlsx"
                with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                    for category, group in df_sorted.groupby('Category'):
                        group.to_excel(writer, sheet_name=str(category)[:31], index=False)
                print(f"Saved {output_file}")
            else:
                print(f"'Category' column not found in {xlsx_file}. Skipping file.")
        else:
            print(f"'Genes' column not found in {xlsx_file}. Skipping file.")

# Example usage:
# directory_path = "/path/to/csv/files"
# lookup_dict = {"FBgn000001": "geneA", "FBgn000002": "geneB"}  # Replace with actual mappings
# process_csv_files(directory_path, lookup_dict)
