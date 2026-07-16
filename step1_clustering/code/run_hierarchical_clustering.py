import os, re, warnings
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings('ignore')

LUIS = '/home/kxj190026/scratch/others/camelomics/luis_data/Donors_included_after_biopsy_QCed'
RENAL = '/home/kxj190026/scratch/others/jeremy/procurement-biopsy-pathomics-ml/data/Renal_Data.csv'
OUTDIR = '/home/kxj190026/scratch/others/jeremy/procurement-biopsy-pathomics-ml/data'
RS = 382025
K_LOCAL = 4
COLOR_KW = ['Blue','Red','Green','Color','Stain']

TISSUES = {
    'glomeruli': 'final_combined_features_non_globally_sclerotic_gloms.xlsx',
    'tubules': 'final_combined_features_tubules.xlsx',
    'arteries': 'final_combined_features_arteriesarterioles.xlsx',
    'scl_gloms': 'final_combined_features_globally_sclerotic_gloms.xlsx'
}

# Normalize folder/slide names for matching (from old script)
def cid(s):
    s = str(s).strip().replace(' - ',' ').replace('_',' ')
    s = re.sub(r' -(\d+)$', r' \1', s)
    s = re.sub(r' PAS (\d+)', r' \1', s, flags=re.I)
    s = re.sub(r' PAS$', '', s, flags=re.I)
    return re.sub(r' +', ' ', s).strip()

# Build folder lookup
dirs = [d for d in os.listdir(LUIS)
        if not d.endswith('.xlsx') and os.path.isdir(os.path.join(LUIS,d))]
lc = {}
for d in dirs:
    k = cid(d); lc[k] = d
    kn = re.sub(r' \d+$','',k).strip()
    if kn not in lc: lc[kn] = d

def find_folder(slide):
    k = cid(slide)
    if k in lc: return lc[k]
    return lc.get(re.sub(r' \d+$','',k).strip(), None)

renal = pd.read_csv(RENAL)
slide_numbers = renal['Slide_number'].tolist()
print(f"Total patients: {len(slide_numbers)}")

# Test matching
matched = sum(1 for s in slide_numbers if find_folder(s) is not None)
print(f"Slides matched to folders: {matched}/{len(slide_numbers)}")

def get_patient_data(tissue_file):
    patient_data = {}
    for slide in slide_numbers:
        folder = find_folder(slide)
        if folder is None:
            continue
        fp = os.path.join(LUIS, folder, tissue_file)
        if not os.path.exists(fp):
            continue
        try:
            df = pd.read_excel(fp)
            if 'compartment_id' in df.columns:
                df = df.drop(columns=['compartment_id'])
            numeric = df.select_dtypes(include=[np.number])
            keep = [c for c in numeric.columns
                    if not any(k in c for k in COLOR_KW)]
            numeric = numeric[keep]
            if len(numeric) > 0:
                patient_data[slide] = numeric.values
        except:
            pass
    return patient_data

def hierarchical_cluster_averaging(patient_data, k_local, tissue_name):
    sample = next(iter(patient_data.values()))
    n_features = sample.shape[1]
    results = {}
    zero_count = 0
    total_entries = 0

    for slide, objects in patient_data.items():
        objects = np.nan_to_num(objects, nan=0.0)
        n_objects = len(objects)
        row = np.zeros(k_local * n_features)

        if n_objects == 0:
            zero_count += k_local * n_features
        elif n_objects < k_local:
            k_use = n_objects
            objects = np.nan_to_num(objects, nan=0.0)
            km = KMeans(n_clusters=k_use, random_state=RS, n_init=5)
            labels = km.fit_predict(objects)
            for c in range(k_use):
                mask = labels == c
                if mask.any():
                    row[c*n_features:(c+1)*n_features] = objects[mask].mean(0)
            zero_count += (k_local - k_use) * n_features
        else:
            scaler = StandardScaler()
            objects_scaled = scaler.fit_transform(objects)
            km = KMeans(n_clusters=k_local, random_state=RS, n_init=10)
            labels = km.fit_predict(objects_scaled)
            for c in range(k_local):
                mask = labels == c
                if mask.any():
                    row[c*n_features:(c+1)*n_features] = objects[mask].mean(0)
                else:
                    zero_count += n_features

        total_entries += k_local * n_features
        results[slide] = row

    zero_pct = 100 * zero_count / total_entries if total_entries > 0 else 0
    print(f"  {tissue_name} k={k_local}: {len(results)} patients, "
          f"zero imputation: {zero_pct:.1f}%")

    col_names = [f'{tissue_name}_hc{c+1}_feat{f+1}'
                 for c in range(k_local) for f in range(n_features)]
    df = pd.DataFrame.from_dict(results, orient='index', columns=col_names)
    df.index.name = 'Slide_number'
    return df.reset_index()

# Process all tissues
all_dfs = []
for tissue_name, tissue_file in TISSUES.items():
    print(f"\nProcessing {tissue_name}...")
    patient_data = get_patient_data(tissue_file)
    print(f"  Patients with data: {len(patient_data)}")
    if len(patient_data) == 0:
        print(f"  Skipping {tissue_name} - no data found")
        continue
    df = hierarchical_cluster_averaging(patient_data, K_LOCAL, tissue_name)
    all_dfs.append(df)

# Merge all tissues with renal slide numbers
print("\nMerging all tissues...")
renal_slides = renal[['Slide_number']].copy()
merged = renal_slides.copy()
for df in all_dfs:
    merged = merged.merge(df, on='Slide_number', how='left')
merged = merged.fillna(0)

print(f"Final shape: {merged.shape}")
non_zero = (merged.iloc[:,1:] != 0).sum().sum()
total = merged.shape[0] * (merged.shape[1]-1)
zero_rate = 100 * (1 - non_zero/total)
print(f"Overall zero imputation: {zero_rate:.1f}%")
print(f"Non-zero values: {non_zero}/{total}")

outfile = os.path.join(OUTDIR, f'Renal_Data_hierarchical_k{K_LOCAL}.csv')
merged.to_csv(outfile, index=False)
print(f"Saved: {outfile}")
