import pandas as pd




#for Neuro experiments:
df_targets = pd.read_csv("/pollard/data/projects/sdrusinsky/pollard_lab/data/enformer_df_targets.csv")
cols = df_targets['description'].tolist()


relevant_cols = []
search_strings = ["brain","neur","glia","cortex","cerebellum","cerebrum","cerebral","parietal lobe","amygdala","spinal","astrocyte",'temporal lobe','frontal lobe','occipital lobe']
for col in cols:
    for search_string in search_strings:
        if search_string.lower() in col.lower():
            if "renal" not in col:
                relevant_cols.append(col)
relevant_cols = [x for x in relevant_cols if x!= ''] #remove emptry strings
relevant_cols = [x.lower() for x in relevant_cols] #coerce all strings to be lowercase
with open("./NeuroEnformerCols.txt",'w') as f:
	for col in relevant_cols:
		f.write(f"{col}\n")
# with open("./NeuroEnformerCols.txt",'r') as f:
# 	file = f.read()
# 	relevant_cols = file.split('\n')




#for Heart experiments:
df_targets = pd.read_csv("/pollard/data/projects/sdrusinsky/pollard_lab/data/enformer_df_targets.csv")
cols = df_targets['description'].tolist()


relevant_cols = []
search_strings = ["heart","atrium","atria","ventricle","blood","cardio",'cardiomyocyte','mural','pericyte','endothelial','lymphoid','myeloid','carotid','mesothelial','vascular','vasculature','artery','vein','vena cava','coronary','stromal','smooth muscle'] #many of these key words came from fig 1 of https://www.nature.com/articles/s41586-020-2797-4#Fig1
for col in cols:
    for search_string in search_strings:
        if search_string.lower() in col.lower():
                relevant_cols.append(col)
relevant_cols = [x for x in relevant_cols if x!= ''] #remove emptry strings
relevant_cols = [x.lower() for x in relevant_cols] #coerce all strings to be lowercase


with open("./HeartEnformerCols.txt",'w') as f:
	for col in relevant_cols:
		f.write(f"{col}\n")
# with open("./HeartEnformerCols.txt",'r') as f:
# 	file = f.read()
# 	relevant_cols = file.split('\n')
# relevant_cols = [x for x in relevant_cols if x!= ''] #remove emptry strings