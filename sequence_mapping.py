from multiprocessing.pool import Pool
from multiprocessing import Process
import pandas as pd
import os
from collections import Counter
import sys

def clean_data(fastq_1_path,fastq_2_path,output_dir,sample_index):
    cleaned_fastq_1 = f'{output_dir}/{sample_index}_filtered_P_1.fq.gz'
    cleaned_fastq_2 = f'{output_dir}/{sample_index}_filtered_P_2.fq.gz'
    if os.path.exists(cleaned_fastq_1) and os.path.exists(cleaned_fastq_1):
        print(f'{sample_index} already cleaned, starting mapping.......')
    else:
        print(f'starting clean sample {sample_index}')
        return_code =os.system(f'trimmomatic PE {fastq_1_path} {fastq_2_path} {output_dir}/{sample_index}_filtered_P_1.fq.gz {output_dir}/{sample_index}_filtered_U_1.fq.gz\
        {output_dir}/{sample_index}_filtered_P_2.fq.gz {output_dir}/{sample_index}_filtered_U_2.fq.gz TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:51')
        if return_code == 0:
            print(f'mession {sample_index} clean finished')
        else:
            print(f'mession {sample_index} clean failed') 
            return 0
    return cleaned_fastq_1,cleaned_fastq_2

def map_star(fstq_1_path,fstq_2_path,out_dir,sample_name):
    out_bam_path = f'{out_dir}Aligned.toTranscriptome.out.bam'
    if os.path.exists(out_bam_path):
        print(f'{out_bam_path} alreading mapped, starting rsem counting......')
    else:
        print(f'starting generating mapping files to {star_our_dir}')
        return_code = os.system(f'STAR --runThreadN 10 --genomeDir "/home/home1/lTang/GRCm39_110_star/" --readFilesIn {fstq_1_path} {fstq_2_path} --readFilesCommand gunzip -c --quantMode TranscriptomeSAM  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {out_dir}')
        if return_code == 0:
            print(f'mession {sample_name} mapping finished')
        else:
            print(f'mession {sample_name} mapping failed') 
            return 0
    return out_bam_path 

def count_rsem(bam_path,resm_ref_dir,output_dir,sample_name):
    if os.path.exists(f'{output_dir}/{sample_name}.genes.results'):
        print(f'{sample_name} already counted, go to the next sample......')
        return f'{output_dir}/{sample_name}.genes.results'
    else:
        print(f'starting rsem counting {sample_name}')
        return_code = os.system(f'rsem-calculate-expression -p 20 --paired-end --bam --append-names --alignments {bam_path} {resm_ref_dir} {output_dir}/{sample_name}')
        if return_code == 0:
            print(f'mession {sample_name} rsem counts finished')
            return f'{output_dir}/{sample_name}.genes.results'
        else:
            print(f'mession {sample_name} rsem counts failed')
            return 0 

def get_expression_matrix(rsem_path,sample_name,matrix_out_dir):
    expresion_out_path = f'{matrix_out_dir}/{sample_name}_expression_counts.csv'
    matrix = pd.read_table(rsem_path,sep='\t')
    matrix.index = matrix['gene_id'].str.split('_').str[1].tolist()
    matrix = matrix.loc[:,['expected_count','TPM','FPKM']]
    matrix.to_csv(expresion_out_path)
    pass

all_sample = []
all_files_path = []
for root,dir,files in os.walk('/home/home1/lTang/BigSpace/mouse_v1_raw/mouse_v1_bulk/'):
    all_sample = all_sample+files
    for name in files:
        all_files_path.append(os.path.join(root,name))
all_sample = [i.split('.')[0][0:-2] for i in all_sample]
for i in Counter(all_sample).keys():
    if Counter(all_sample).get(i) !=2:
        print(f'sample {i} have doublets, checking it and try again')
        sys.exit()
# create "sample : file_path dict"
unique_sample = list(set(all_sample))
sample_path_dict = dict()
for i in unique_sample:
    temp_fast_1 = []
    temp_fast_2 = []
    for j in all_files_path:
        fastq_name = os.path.split(j)[1].split('.')[0]
        if fastq_name == i+'_1':
            temp_fast_1.append(j)
        if fastq_name == i+'_2':
            temp_fast_2.append(j)
    if (len(temp_fast_1)!=1) or (len(temp_fast_1)!=1):
        print(f'sample {i} is duplicated or not found pleace check the path {temp_fast_1} or {temp_fast_1}')
        sys.exit() 
    sample_path_dict[i] = [temp_fast_1[0],temp_fast_2[0]]
# clean each sample and mapping
clean_data_out_dir = './clean_data'
star_out_root = './star_outs'
rsem_out_dir = './rsem_outs'
matrix_out_dir = './matrix_outs'
for path_index in [clean_data_out_dir,star_out_root,rsem_out_dir,matrix_out_dir]:
    if os.path.exists(path_index):
        print(f'directory {path_index} already exist, please check and store it to avoid overwriting')
    else:
        os.mkdir(path_index)

for sample_index in sample_path_dict.keys():
    clean_code= clean_data(sample_path_dict.get(sample_index)[0],sample_path_dict.get(sample_index)[1],clean_data_out_dir,sample_index)
    if clean_code == 0:
        continue
    else:
        temp_clean_fastq_1,temp_clean_fastq_2 = clean_code[0],clean_code[1]
    star_our_dir = os.path.join(star_out_root,sample_index)
    map_return_code = map_star(temp_clean_fastq_1,temp_clean_fastq_2,star_our_dir,sample_index)
    if map_return_code == 0:
        continue
    else:
        temp_bam_path = map_return_code
    rsem_code = count_rsem(temp_bam_path,"/home/home1/lTang/GRCm39_110_rsem/mouse_rsef_reference",rsem_out_dir,sample_index)
    if rsem_code == 0:
        print(f'rsem countting in sample {sample_index} meet some problem, jump to next sample')
        continue
    get_expression_matrix(rsem_code,sample_index,matrix_out_dir)
