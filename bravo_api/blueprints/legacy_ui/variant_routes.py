"""@package Bravo Variant Routes
Provide convenient endpoints for variant queries.
Primary responsibilities are:
    - Consuming arguments from UI
    - Providing routes that use view arguments
    - Wrapping data in web responses.
"""
from flask import Blueprint, make_response, jsonify, send_file, abort, request
from webargs import fields
from marshmallow import validate
from bravo_api.blueprints.legacy_ui import pretty_api, common
import re
import requests

from bs4 import BeautifulSoup
import math
import json

bp = Blueprint('variant_routes', __name__)

parser = common.Parser()


@bp.route('/qc/api')
def qc():
    result = pretty_api.get_qc()

    response = make_response(jsonify(result), 200)
    response.mimetype = 'application/json'
    return(response)


variant_argmap = {
    'variant_id': fields.Str(required=True, validate=validate.Length(min=1),
                             error_messages=common.ERR_EMPTY_MSG)
}


#HX
@bp.route('/tempdata/<string:variant_id>', methods=['GET'])
def get_data(variant_id):
    url = f'https://pheweb.org/UKB-TOPMed/variant/{variant_id}'
    response = requests.get(url)
    return response.text

@bp.route('/tempdata_mod/UKB-TOPMed/<string:variant_id>', methods=['GET'])
def get_data_mod_ukb(variant_id):
    url = f'https://pheweb.org/UKB-TOPMed/variant/{variant_id}'
    response = requests.get(url)
    var_name = "window.variant"
    raw_data = ""
    soup = BeautifulSoup(response.text, 'html.parser')
    scripts = soup.find_all('script', {'type': 'text/javascript'})
    for script in scripts:
        script_code = script.string
        if (script_code) and (var_name in script_code):
            match = re.search(f'{var_name} = (.*?\\}});', script_code)
            if match:
                raw_data = match.group(1)
    if raw_data:
        data_dict = json.loads(raw_data)
        result = {
            "data": [],
            "lastPage": None,
            "meta": {"build": ["GRCh38"]}
        }
        record_id = 1

        for pheno in data_dict["phenos"]:
            entry = {
                "beta": pheno["beta"],
                "chromosome": data_dict["chrom"],
                "position": data_dict["pos"],
                "id": record_id,
                "ref_allele": data_dict["ref"],
                "variant": data_dict["variant_name"],
                "log_pvalue": -math.log10(pheno["pval"]),
                "trait_group": pheno["category"],
                "trait": pheno["phenocode"],
                "trait_label": pheno["phenostring"],
                "num_cases": pheno["num_cases"],
            }
            result["data"].append(entry)
            record_id += 1
        json_result = json.dumps(result)
        return json_result
    else:
        return json.dumps({})

@bp.route('/tempdata_mod/freeze5/<string:variant_id>', methods=['GET'])
def get_data_mod_f5(variant_id):
    url = f'http://r5.finngen.fi/variant/{variant_id}'
    response = requests.get(url)
    var_name = "window.results"
    var_name2 = "window.variant"
    soup = BeautifulSoup(response.text, 'html.parser')
    scripts = soup.find_all('script', {'type': 'text/javascript'})
    raw_data2 = ""
    raw_data_var = ""
    for script in scripts:
        script_code = script.string
        if (script_code) and (var_name in script_code):
            match = re.search(f'{var_name} = (.*?\\]);', script_code)
            if match:
                raw_data2 = match.group(1)
        if (script_code) and (var_name2 in script_code):
            # print("yes")
            match = re.search(f'{var_name2} = (.*?\\}});', script_code)
            if match:
                raw_data_var = match.group(1)
    if raw_data2 and raw_data_var:
        data_dict = json.loads(raw_data2)
        data_var = json.loads(raw_data_var)
        
        # if not isinstance(data_dict, list):
        #     return json.dumps({})
        # if not isinstance(data_var, dict):
        #     return json.dumps({})
        
        result = {
            "data": [],
            "lastPage": None,
            "meta": {"build": ["GRCh38"]}
        }
        record_id = 1

        for pheno in data_dict:
            try:
                pval_float = float(pheno["pval"])
            except (ValueError, TypeError):
                pval_float = None
            
            entry = {
                "beta": pheno["beta"],
                "chromosome": data_var["chr"],
                "position": data_var["pos"],
                "id": record_id,
                "ref_allele": data_var["ref"],
                "variant": data_var["varid"],
                "log_pvalue": -math.log10(pheno["pval"]) if pval_float is not None else None,
                "trait_group": pheno["category"],
                "trait": pheno["phenocode"],
                "trait_label": pheno["phenostring"]
            }
            result["data"].append(entry)
            record_id += 1
        json_result = json.dumps(result)
        return json_result
    else:
        return json.dumps({})
        

@bp.route('/test', methods=['GET'])
def test_route():
    return "Test is working!"
# HX

@bp.route('/variant/api/snv/<string:variant_id>')
@parser.use_kwargs(variant_argmap, location='view_args')
def variant(variant_id):
    result = pretty_api.get_variant(variant_id)

    response = make_response(jsonify(result), 200)
    response.mimetype = 'application/json'
    return(response)


@bp.route('/variant/api/snv/cram/summary/<string:variant_id>')
@parser.use_kwargs(variant_argmap, location='view_args')
def variant_cram_info(variant_id):
    result = pretty_api.get_variant_cram_info(variant_id)

    response = make_response(jsonify(result), 200)
    response.mimetype = 'application/json'
    return(response)


variant_cram_argmap = {
    'variant_id': fields.Str(required=True, validate=validate.Length(min=1),
                             error_messages=common.ERR_EMPTY_MSG),
    'sample_no': fields.Int(required=True, validate=validate.Range(min=1),
                            error_messages=common.ERR_GT_ZERO_MSG),
    'sample_het': fields.Bool(required=True)
}


@bp.route('/variant/api/snv/cram/<string:variant_id>-<int:sample_het>-<int:sample_no>')
@parser.use_kwargs(variant_cram_argmap, location='view_args')
def variant_cram(variant_id, sample_het, sample_no):
    range_header = request.headers.get('Range', None)
    start = None
    stop = None
    if range_header:
        m = re.search(r'(\d+)-(\d*)', range_header)
        if m:
            start = int(m.group(1))
            stop = int(m.group(2))

    result = pretty_api.get_variant_cram(variant_id, sample_het, sample_no, start, stop)
    if result is None:
        print(f'start: {start} stop: {stop}')
        abort(404)

    response = make_response(result['file_bytes'], 206)
    response.headers['Content-Range'] = f'bytes {result["start_byte"]}-{result["stop_byte"]}/{result["file_size"]}'
    response.mimetype = 'application/octet-stream'
    response.direct_passthrough = True

    return(response)


@bp.route('/variant/api/snv/crai/<string:variant_id>-<int:sample_het>-<int:sample_no>')
@parser.use_kwargs(variant_cram_argmap, location='view_args')
def variant_crai(variant_id, sample_no, sample_het):
    result = pretty_api.get_variant_crai(variant_id, sample_no, sample_het)
    if result is None:
        abort(404)
    response = make_response(send_file(result, as_attachment=False))
    return response
