import pandas as pd
import numpy as np
from pandas import HDFStore
import h5py
import os
from pathlib import Path
import shutil

oneKG_vars = "/pollard/data/projects/sdrusinsky/pollard_lab/data/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf"
TMPDIR = os.environ.get("TMPDIR") #get the location of the local scratch folder on the compute node. This is where output will be saved to and then moved to pollard server after the fact
scratch_path = f"{TMPDIR}/get_rare_1kg_vars/" #where temporary will be saved
Path(scratch_path).mkdir(parents=True, exist_ok=True) #create folder structure in local scratch that data will be saved in
assert os.path.exists(scratch_path)

trait='Rare1KGVariants'
global_output_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants"
def save_vars_as_vcf(study_tag_variants,output_dir,trait):

    
    variant_df = pd.DataFrame(columns = ["CHROM","POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
    variant_df = variant_df.append(study_tag_variants.rename(columns = {'#CHROM':'CHROM'})[['CHROM','POS','REF','ALT','INFO']])

    variant_df = variant_df.drop_duplicates(['CHROM','POS','REF','ALT'])

    #fill empty columns
    variant_df['ID'] = ""
    variant_df['QUAL'] = ""
    variant_df['FILTER'] = ""

    #force column order
    variant_df = variant_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']] 
    
    #save vcf with no header, since Sei's parser does not seem to like it
    sei_VCF = os.path.join(output_dir,f"PsychENCODE_GWASVariants_{trait}_SeiNoHeader.vcf")
    variant_df.to_csv(sei_VCF, sep="\t",index=False,header=False)

    #save VCF with header for Enformer and Akita
    header = """##fileformat=VCFv4.1
    ##fileDate=20220819
    ##source=/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/DNVs/get_rare_1kg_vars.py
    ##reference=/wynton/home/hernandez/shirondru/pollard_lab/data/hg38_genome.fa
    #CHROM POS ID REF ALT QUAL FILTER INFO
    """
    #write VCF header
    #save chunks of df to later parallelize model predictions on each chunk
    output_VCF = os.path.join(output_dir,f"PsychENCODE_GWASVariants_{trait}.vcf")
    with open(output_VCF, 'w') as vcf:
        vcf.write(header)

    variant_df.to_csv(output_VCF, sep="\t",mode = 'a',index=False) #append data below vcf header


colnames = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]
idx = 0
chunksize = 100000
for chunk in pd.read_csv(oneKG_vars,skiprows = 41, names = colnames,chunksize = chunksize,sep='\t'):
    chunk = chunk[chunk['FILTER']=='PASS']
    
    #get population allele freq data in the INFO column
    data = chunk['INFO'].str.split(';',expand=True)[[3,4,5,6,7,8]]
    
    
    # remove string characters and keep only float values
    for population in range(3,9):
        data[population] = data[population].str.extract('(0\.\d+|0)').astype(float)
    #add allele frequency float columns back to the chunk containing variant position and allele nucleotides
    var_data = pd.concat([chunk.drop(columns = ['FILTER','QUAL']),data],axis=1)
    
    ### Keep only variants with MAF >=0.05 in any population. Only keep columns related to variant position and allele nucleotides
    ### append to this dataframe while iterating through chunks
    if idx == 0:
        vars_to_keep = var_data[(var_data[[3,4,5,6,7,8]] <= 0.01).all(1)][['#CHROM','POS','ID','REF','ALT','INFO']]
        idx += 1
    else:
        vars_to_keep = vars_to_keep.append(var_data[(var_data[[3,4,5,6,7,8]] <= 0.01).all(1)][['#CHROM','POS','ID','REF','ALT','INFO']])
        idx+=1
    if idx % 100 == 0:
        print(f"{(chunksize*idx)/73257674}% finished")

random_vars_to_keep = vars_to_keep.sample(1000000)
save_vars_as_vcf(random_vars_to_keep,scratch_path,trait)
result = os.path.exists(f"{scratch_path}/PsychENCODE_GWASVariants_{trait}_SeiNoHeader.vcf")
print(f"PsychENCODE_GWASVariants_{trait}_SeiNoHeader.vcf was saved properly in scratch: {result}")
result = os.path.exists(f"{scratch_path}/PsychENCODE_GWASVariants_{trait}.vcf")
print(f"PsychENCODE_GWASVariants_{trait}.vcf was saved properly in scratch: {result}")

shutil.move(f"{scratch_path}/PsychENCODE_GWASVariants_{trait}_SeiNoHeader.vcf",f"{global_output_dir}/PsychENCODE_GWASVariants_{trait}_SeiNoHeader.vcf")
shutil.move(f"{scratch_path}/PsychENCODE_GWASVariants_{trait}.vcf",f"{global_output_dir}/PsychENCODE_GWASVariants_{trait}.vcf")

result = os.path.exists(f"{global_output_dir}/PsychENCODE_GWASVariants_{trait}_SeiNoHeader.vcf")
print(f"PsychENCODE_GWASVariants_{trait}_SeiNoHeader.vcf was saved properly in global outdir: {result}")
result = os.path.exists(f"{global_output_dir}/PsychENCODE_GWASVariants_{trait}.vcf")
print(f"PsychENCODE_GWASVariants_{trait}.vcf was saved properly in global outdir: {result}")

