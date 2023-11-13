import pandas as pd
import os

def load_and_update_tsv(file_path):
    """
    Function to load a tsv file, perform replacements based on predefined dictionaries,
    and write the updated dataframe back to the file
    """
    
    # Check if the file exists
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return

    # Load the tsv file
    df = pd.read_csv(file_path, sep='\t')

    # Define a dictionary for the replacements, organized alphabetically by body part
    replacement_dict_curate = {
    
    # Abdomen replacements
    "apis mellifera rectum": "rectum",
    "ovary": "ovaries",
    "eviscerated abdomen with attached fat bodies": "eviscerated abdomen and fat body",
    "eviscerated abdomen with attached fat bodies and active ovaries": "eviscerated abdomen with attached fat bodies and ovaries",
    
    # Antennae replacements
    "4 pairs of antennae": "antennae",
    "antenna": "antennae",
    
    # Brain replacements
    "whole brain": "brain",
    "5 pooled brains": "brain",
    "brain (pool of 4 brains)": "brain",
    "optic lobes (ol) region of the brain": "optic lobes",
    
    # Reproductive replacements
    "drone germline": "germline",
    "queen germline": "germline",
    "drone thorax": "thorax",
    "spermetheca_virgin": "spermetheca",
    "spermethecia_mated": "spermetheca",
    "testicles": "testis",
    
    # Embryo replacements
    "drone embryo": "embryo",
    "embryos": "embryo",
    "whole embryos": "embryo",
    "worker embryo": "embryo",
    
    # Gland replacements
    "malpighian tubule": "malpighian tubules",
    "mandibular gland": "mandibular glands",
    "Mandibular_Gland_Apis_mellifera_Nurse_C1": "mandibular glands",
    "Mandibular_Gland_Apis_mellifera_Nurse_C3": "mandibular glands",
    "Mandibular_Gland_Nurse_C2": "mandibular glands",
    "hypopharyngeal glands tissue": "hypopharyngeal glands",
    "hypopharyngeal gland": "hypopharyngeal glands",
    "venom system": "venom gland",
    
    # Miscellaneous replacements
    "whole abdomen": "abdomen",
    "fatbody": "fat body",
    "gonadal": "male gonad",
    "gut tissue": "gut",
    "heads": "head",
    "pupa": "pupae",
    "rna of the whole pupae": "pupae",
    "intestine": "gut",
    "pooled pylori from 7 honey bees": "pylori",
    "membranous wings": "wings",
    "second thoracic ganglia": "2nd thoracic ganglia",
    "tissue": "unknown",
    "not applicable": "unknown",
    "unknow": "unknown",
    "": "unknown",
    
    # Larvae
    "larvea": "larvae",
    "whole body larvae": "larvae",    
    "larva": "larvae",
    "larvae body": "larvae",
    "whole larva": "larvae",
    "whole larvae": "larvae",
    "whole larva body": "larvae",
    "whole larvae inoculated with b8 rep 1": "larvae",    
    "whole larvae inoculated with b8 rep 2": "larvae",
    "whole larvae inoculated with b8 rep 3": "larvae",
    "whole larvae inoculated with c6 rep 1": "larvae",
    "whole larvae inoculated with c6 rep 2": "larvae",    
    "whole larvae inoculated with c6 rep 3": "larvae",
    
    # Mushroom body replacements
    "mushroom bodies": "mushroom body",
    "mushroom bodies (mb) region of the brain": "mushroom body",
    "mushroom bodies of brain": "mushroom body",
    "mushroom body (dissected from brain)": "mushroom body",
    
    # Whole body replacements
    "whole bodies": "whole adult body",
    "whole body": "whole adult body",
    "whole body homogenates": "whole adult body",
    "10 pooled whole bodies": "whole adult body",
    "body": "whole adult body",
    "entirety": "whole adult body",
    "total body": "whole adult body",
    "total individual": "whole adult body",    
    "total individuals": "whole adult body",
    "total body": "whole adult body",
    "total body": "whole adult body",
    "whole": "whole adult body",
    "whole adult female": "whole adult body",
    "whole bee": "whole adult body",
    "whole insect": "whole adult body",
    "whole tissue": "whole adult body",
    "whole worker bees": "whole adult body",
    "whole_body_1": "whole adult body",    
    "whole_body_2": "whole adult body",    
    "whole_body_3": "whole adult body",
    "whole_body_4": "whole adult body",
    "whole_body_5": "whole adult body",
    "whole_body_6": "whole adult body",
    "whole_body_7": "whole adult body",
    "whole_body_8": "whole adult body",
    "whole_body_9": "whole adult body",
    "whole_body_10": "whole adult body",    
    "whole_body_11": "whole adult body", 
    "whole_body_12": "whole adult body",    
    "whole_body_13": "whole adult body",
    "whole_body_14": "whole adult body",
    "whole_body_15": "whole adult body",
    "whole_body_16": "whole adult body",
    "whole_body_17": "whole adult body",
    "whole_body_18": "whole adult body",
    "whole_body_19": "whole adult body",  
    "whole_body_20": "whole adult body",    
    "whole_body_21": "whole adult body",    
    "whole_body_22": "whole adult body",    
    "whole_body_23": "whole adult body",
    "whole_body_24": "whole adult body",
    "whole_body_25": "whole adult body",
    "whole_body_26": "whole adult body",
    "whole_body_27": "whole adult body",
    "whole_body_28": "whole adult body",
    }

    replacement_dict_scientific_name = {
        "": "Apis mellifera",
        "Apis mellifera capensis": "Apis mellifera",
        "Apis mellifera carnica": "Apis mellifera",
        "Apis mellifera intermissa": "Apis mellifera",
        "Apis mellifera ligustica": "Apis mellifera",
        "Apis mellifera mellifera": "Apis mellifera",
        "Apis mellifera scutellata": "Apis mellifera",
        "Apis mellifera syriaca": "Apis mellifera",
        "Apis mellifrea": "Apis mellifera",
        }

        # Replace the values in the "curate_group" column
    df['curate_group'] = df['curate_group'].replace(replacement_dict_curate)
    df['scientific_name'] = df['scientific_name'].replace(replacement_dict_scientific_name)

        # Write the data back to the .tsv file
    try:
        df.to_csv(file_path, sep='\t', index=False)
    except Exception as e:
        print(f"Error writing to file: {e}")

def main():
    file_path = './metadata/metadata.tsv'
    load_and_update_tsv(file_path)

if __name__ == "__main__":
    main()

