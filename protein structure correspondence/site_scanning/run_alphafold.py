import numpy as np
import os
import shutil
import multiprocessing as mp



path = 'path_to_fasta'
paths=[]
for file in os.listdir(path):
    filename=os.fsdecode(file)
    if filename.endswith('.pdb'):
        paths.append(f'{path}/{filename}')


def alphafold(fa,protein_name=1):

   


    DATA='/rds/project/rds-dMMtPvqHWv4/data' 
    '''run monomer alphafold'''
#     os.system(f'run_alphafold \
# 		--bfd_database_path={DATA}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
#         --pdb70_database_path=/data/pdb70/pdb70 \
# 		--uniclust30_database_path={DATA}/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
# 		--output_dir={...} \
# 		--fasta_paths={fa} \
# 		--max_template_date=2020-05-14  \
# 		--db_preset=full_dbs \
# 		--use_gpu_relax=True \
# 		--cpus 32')
    '''run multimer alphafold'''
    SINGULARITY_DIR='/usr/local/Cluster-Apps/singularity/images/'
    os.system(f'singularity run --env TF_FORCE_UNIFIED_MEMORY=1,XLA_PYTHON_CLIENT_MEM_FRACTION=4.0,OPENMM_CPU_THREADS=32 -B $DATA:{DATA} -B .:/etc --pwd /app/alphafold \
              --nv {SINGULARITY_DIR}/alphafold-2.1.2.sif \
              --benchmark   \
            --output_dir {...} \
            --fasta_paths {fa}    \
            --data_dir {DATA} \
            --max_template_date=2020-05-14   \
            --uniref90_database_path {DATA}/uniref90/uniref90.fasta   \
            --mgnify_database_path {DATA}/mgnify/mgy_clusters_2018_12.fa   \
            --uniclust30_database_path {DATA}/uniclust30/uniclust30_2018_08/uniclust30_2018_08   \
            --bfd_database_path {DATA}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt   \
            --template_mmcif_dir {DATA}/pdb_mmcif/mmcif_files   \
            --obsolete_pdbs_path {DATA}/pdb_mmcif/obsolete.dat   \
            --db_preset full_dbs   \
            --model_preset multimer \
            --pdb_seqres_database_path {DATA}/pdb_seqres/pdb_seqres.txt \
            --uniprot_database_path {DATA}/uniprot/uniprot.fasta \
            --cpus 32')
        
    


pool = mp.Pool(4)
pg = [pool.apply_async(alphafold,args=(i,1))for i in paths]
pg = [p.get() for p in pg]
