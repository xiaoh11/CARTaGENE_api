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

# def fetch_vcf_positions(name):
#     # print("fetch is called")
#     positions = []
#     vcf_file = current_app.config['CLINVAR_VCF']
    
#     gene = variants.get_gene(name, True)
#     if gene is None:
#         return positions
#     chromosome = gene['chrom']
#     start = gene['start']
#     end = gene['stop']
    
#     with pysam.VariantFile(vcf_file) as vcf:
#         for record in vcf.fetch(chromosome, start, end):
#             base_info = f"{record.chrom}-{record.pos}-{record.ref}-"
#             clnsig = record.info.get('CLNSIG')
#             if record.alts:
#                 for alt in record.alts:
#                     variant_info = base_info + str(alt)
#                     # positions.append(variant_info)
#             else:
#                 positions.append([base_info, clnsig]) 
#             if clnsig:
#                 for sig in clnsig:
#                     positions.append([variant_info, record.id, str(sig)])
    
#     return positions
def fetch_vcf_positions(name):
    positions = []
    vcf_file = current_app.config['CLINVAR_VCF']
    
    gene = variants.get_gene(name, True)
    if gene is None:
        return positions

    chromosome = gene['chrom']
    start = gene['start']
    end = gene['stop']

    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch(chromosome, start, end):
            base_info = f"{record.chrom}-{record.pos}-{record.ref}-"
            clnsig = record.info.get('CLNSIG', [])

            variant_infos = []   
            if record.alts:
                for alt in record.alts:
                    variant_info = base_info + str(alt)
                    variant_infos.append(variant_info) 
            else:
                variant_info = base_info + "?"  #   
                variant_infos.append(variant_info)

            for variant_info in variant_infos:
                if clnsig:
                    for sig in clnsig:
                        positions.append([variant_info, record.id, str(sig)])
                else:
                    positions.append([variant_info, "NA", "Unknown"])  # 

    return positions


def fetch_vcf_positions_by_region(chrom, start, stop):
    print("fetch is called")
    positions = []
    vcf_file = current_app.config['CLINVAR_VCF']
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch(chrom, start, stop):
            base_info = f"{record.chrom}-{record.pos}-{record.ref}-"
            clnsig = record.info.get('CLNSIG', [])

            variant_infos = []   
            if record.alts:
                for alt in record.alts:
                    variant_info = base_info + str(alt)
                    variant_infos.append(variant_info) 
            else:
                variant_info = base_info + "?"  #   
                variant_infos.append(variant_info)

            for variant_info in variant_infos:
                if clnsig:
                    for sig in clnsig:
                        positions.append([variant_info, record.id, str(sig)])
                else:
                    positions.append([variant_info, "NA", "Unknown"])
    
    return positions


def compare_db_with_clinVar(gene_name):
    gene = variants.get_gene(gene_name, True)
    if gene is None:
        return []

    gene_id = gene['gene_id']
    xpos_start = make_xpos(gene['chrom'], gene['start'])
    xpos_end = make_xpos(gene['chrom'], gene['stop'])

    query = {
        'xpos': {'$gte': xpos_start, '$lte': xpos_end},
        'annotation.genes.name': gene_id  
    }

    projection = {'variant_id': True, '_id': True}  
    cursor = mongo.db.snv.find(query, projection)

    variant_ids = [doc['variant_id'] for doc in cursor]
    
    clinvar_variants = fetch_vcf_positions(gene_name)
    
    for variant_info in clinvar_variants:
        if variant_info[0] in variant_ids:
            variant_info.append(1)
            # overlapped_variants.append(variant_info[0])
        else:
            variant_info.append(0)

    return clinvar_variants

def compare_db_with_clinVar_region(chrom, start, stop):
    xpos_start = make_xpos(chrom, start)
    xpos_end = make_xpos(chrom, stop)

    query = {
        'xpos': {'$gte': xpos_start, '$lte': xpos_end}
    }

    projection = {'variant_id': True, '_id': False}  
    cursor = mongo.db.snv.find(query, projection)

    variant_ids = [doc['variant_id'] for doc in cursor]
    
    clinvar_variants = fetch_vcf_positions_by_region(chrom, start, stop)
    
    for variant_info in clinvar_variants:
        if variant_info[0] in variant_ids:
            variant_info.append(1)
        else:
            variant_info.append(0)
    return clinvar_variants