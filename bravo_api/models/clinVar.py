from bravo_api.models.database import mongo
from bravo_api.models.utils import make_xpos
from flask import current_app
import pymongo
from bson.objectid import ObjectId
from bson.regex import Regex
import functools
from intervaltree import Interval, IntervalTree
from collections import Counter
from bravo_api.models import variants, qc_metrics, sequences
import pysam
import os

# vcf_file_path = "~/scratch/clinvar_20231007.vcf.gz"
def fetch_vcf_positions(name):
    print("fetch is called")
    positions = []
    # vcf_file = "/home/xiaoh11/scratch/clinvar_20231007.vcf.gz"
    # vcf_dir = current_app.config['CLINVAR_DIR']
    # vcf_file = f'{vcf_dir}/clinvar_20231007.vcf.gz'
    vcf_file = current_app.config['CLINVAR_VCF']
    
    gene = variants.get_gene(name, True)
    if gene is None:
        return positions
    print(gene['chrom'], gene['start'], gene['stop'])
    chromosome = gene['chrom']
    start = gene['start']
    end = gene['stop']
    
    with pysam.VariantFile(vcf_file) as vcf:
        # for record in vcf.fetch(chromosome, start, end):
        #     clnsig = record.info.get('CLNSIG')
        #     base_info = f"{record.chrom}-{record.pos}-{record.ref}-"
            
        #     if record.alts:
        #         for alt in record.alts:
        #             variant_info = base_info + str(alt)
        #             positions.append(variant_info)
        #     else:
        #         positions.append([base_info, clnsig]) 
        for record in vcf.fetch(chromosome, start, end):
            base_info = f"{record.chrom}-{record.pos}-{record.ref}-"
            clnsig = record.info.get('CLNSIG')
            if record.alts:
                for alt in record.alts:
                    variant_info = base_info + str(alt)
                    positions.append(variant_info)
            else:
                positions.append([base_info, clnsig]) 
            if clnsig:
                for sig in clnsig:
                    positions.append([variant_info, str(sig)])
    
    return positions

# positions_list = fetch_vcf_positions(vcf_file_path, chromosome, start, end)
# if not os.path.exists(vcf_file_path):
#     print(f"File {vcf_file_path} does not exist!")
#     exit()
# else:
#     print(positions_list)
