from flask import Flask, render_template, request, jsonify, Response, send_file
import pandas as pd
from collections import OrderedDict
import json
import os
import time
import zipfile
import tempfile

app = Flask(__name__, template_folder=".")
DATA_DIR = "."

app.config['JSON_SORT_KEYS'] = False


GR_CAS9 = None
DEEPMENS = None
TRANS = None
GR_ABE = None
GR_CBE = None
DEEPBE_ABE = None
DEEPBE_CBE = None
BEDEE_ABE = None
BEDEE_CBE = None
DEEPBE_ABE_SCORES = None
DEEPBE_CBE_SCORES = None
BEDEE_ABE_SCORES = None
BEDEE_CBE_SCORES = None
CAS9 = None
BE = None


data_loaded = False

def dataframe_to_ordered_records(df):
    """Convert DataFrame to list of ordered dicts, preserving column order"""
    records = []
    for _, row in df.iterrows():
        record = OrderedDict()
        for col in df.columns:
            record[col] = row[col]
        records.append(record)
    return records


data_files = {
    'Prism-CRISPR_cas9': 'Prism-CRISPR_cas9.csv',
    'DeepMEns': 'DeepMEns.csv',
    'TransCrispr': 'TransCrispr.csv',
    'Prism-CRISPR_abe': 'Prism-CRISPR_abe.csv',
    'Prism-CRISPR_cbe': 'Prism-CRISPR_cbe.csv',
    'DeepBE_abe': 'DeepBE_abe.csv',
    'DeepBE_cbe': 'DeepBE_cbe.csv',
    'BEDeepon_abe': 'BEDeepon_abe.csv',
    'BEDeepon_cbe': 'BEDeepon_cbe.csv'
}

# 压缩文件路径
ZIP_FILE = 'data.zip'

def load_csv(filename, filepath):
    """Load CSV file, extracting from zip if needed"""
    # 检查是否存在压缩文件
    if os.path.exists(ZIP_FILE):
        print(f"Extracting {filename} from {ZIP_FILE}")
        with zipfile.ZipFile(ZIP_FILE, 'r') as zip_ref:
            # 检查文件是否在压缩包中
            if filename in zip_ref.namelist():
                # 创建临时目录
                with tempfile.TemporaryDirectory() as tmpdir:
                    # 解压文件到临时目录
                    zip_ref.extract(filename, tmpdir)
                    # 读取解压后的文件
                    extracted_path = os.path.join(tmpdir, filename)
                    df = pd.read_csv(extracted_path)
                    return df
            else:
                print(f"{filename} not found in {ZIP_FILE}")
                # 如果压缩包中没有，尝试直接读取
                if os.path.exists(filepath):
                    print(f"Loading {filename} from local file")
                    df = pd.read_csv(filepath)
                    return df
                else:
                    raise FileNotFoundError(f"{filename} not found in {ZIP_FILE} or locally")
    else:
        # 如果没有压缩文件，直接读取本地文件
        print(f"Loading {filename} from local file")
        df = pd.read_csv(filepath)
        return df

def load_all_data():

    global GR_CAS9, DEEPMENS, TRANS, GR_ABE, GR_CBE, DEEPBE_ABE, DEEPBE_CBE, BEDEE_ABE, BEDEE_CBE
    global DEEPBE_ABE_SCORES, DEEPBE_CBE_SCORES, BEDEE_ABE_SCORES, BEDEE_CBE_SCORES
    global CAS9, BE, data_loaded
    
    if data_loaded:
        return
    
    GR_CAS9 = load_csv('Prism-CRISPR_cas9.csv', data_files['Prism-CRISPR_cas9'])
    DEEPMENS = load_csv('DeepMEns.csv', data_files['DeepMEns'])
    TRANS = load_csv('TransCrispr.csv', data_files['TransCrispr'])
    GR_ABE = load_csv('Prism-CRISPR_abe.csv', data_files['Prism-CRISPR_abe'])
    GR_CBE = load_csv('Prism-CRISPR_cbe.csv', data_files['Prism-CRISPR_cbe'])
    DEEPBE_ABE = load_csv('DeepBE_abe.csv', data_files['DeepBE_abe'])
    DEEPBE_CBE = load_csv('DeepBE_cbe.csv', data_files['DeepBE_cbe'])
    BEDEE_ABE = load_csv('BEDeepon_abe.csv', data_files['BEDeepon_abe'])
    BEDEE_CBE = load_csv('BEDeepon_cbe.csv', data_files['BEDeepon_cbe'])
    

    DEEPBE_ABE_SCORES = prepare_be_auxiliary(DEEPBE_ABE, ["DeepBE_ABE7", "DeepBE_ABE8e", "DeepBE_ABEmax"])
    DEEPBE_CBE_SCORES = prepare_be_auxiliary(DEEPBE_CBE, ["DeepBE_BE4", "DeepBE_CBE4max", "DeepBE_AID"])
    BEDEE_ABE_SCORES = prepare_be_auxiliary(BEDEE_ABE, ["BEDeepon_ABE7", "BEDeepon_ABE8e", "BEDeepon_ABEmax"])
    BEDEE_CBE_SCORES = prepare_be_auxiliary(BEDEE_CBE, ["BEDeepon_BE4", "BEDeepon_CBE4max", "BEDeepon_AID"])
    

    CAS9 = build_cas9()
    BE = build_be()
    
    data_loaded = True

def prepare_be_auxiliary(df, score_cols):
    aux = df.copy()
    aux["target_base"] = aux["sgrna"].astype(str).str[:20].str.upper()
    aux["outcome_base"] = aux["outcome"].astype(str).str[:20].str.upper()
    keep_cols = ["Name", "target_base", "outcome_base", *score_cols]
    return aux[keep_cols].drop_duplicates(subset=["Name", "target_base", "outcome_base"])


def build_cas9():
    base = GR_CAS9.copy()
    base["sgRNA_base"] = base["sgrna"].astype(str).str[:20]

    # only keep score columns from auxiliary datasets to prevent duplicates
    deep_cols = ["sgrna", "DeepMEns_WT", "DeepMEns_ESP", "DeepMEns_HF"]
    trans_cols = ["sgrna", "TransCrispr_WT", "TransCrispr_ESP", "TransCrispr_HF"]

    df = base.merge(DEEPMENS[deep_cols], on="sgrna", how="left")
    df = df.merge(TRANS[trans_cols], on="sgrna", how="left")
    return df


def build_be():
    abe = GR_ABE.copy()
    cbe = GR_CBE.copy()
    df = pd.concat([abe,cbe], ignore_index=True)
    # normalize target base to uppercase 20bp for case-insensitive matching
    df["target_base"] = df["sgrna"].astype(str).str[:20].str.upper()
    df["outcome_base"] = df["outcome"].astype(str).str[:20].str.upper()

    # merge DeepBE and BEDeepon data so that scores are available
    df = df.merge(DEEPBE_ABE_SCORES, on=["Name", "target_base", "outcome_base"], how="left")
    df = df.merge(DEEPBE_CBE_SCORES, on=["Name", "target_base", "outcome_base"], how="left")
    df = df.merge(BEDEE_ABE_SCORES, on=["Name", "target_base", "outcome_base"], how="left")
    df = df.merge(BEDEE_CBE_SCORES, on=["Name", "target_base", "outcome_base"], how="left")
    return df


@app.route("/")
def home():
    return render_template("database.html")


@app.route("/database")
def database():
    return render_template("database.html")


@app.route('/style.css')
def style_css():
    return send_file('style.css')

@app.route('/main.js')
def main_js():
    return send_file('main.js')


@app.route("/api/search", methods=["POST"])
def search():

    load_all_data()
    
    data = request.json
    omim = data.get("omim_id","").strip()
    mutation = data.get("mutation","").strip()
    mode = data.get("mode")
    editor = data.get("editor")
    sort_by = data.get("sort_by")  # optional column to sort results
    
    if mode not in ["cas9","be"] or editor is None:
        return jsonify([])

    # select dataframe
    if mode=="cas9":
        df = CAS9
        if omim:
            # Convert omim to float for comparison, handles both "136880" and "136880.0"
            try:
                omim_val = float(omim)
                df = df[df["OMIM_ID"] == omim_val]
            except:
                df = df[df["OMIM_ID"].astype(str)==omim]
        if mutation:
            # case-insensitive substring match on Name
            df = df[df["Name"].str.contains(mutation, case=False, na=False, regex=False)]
        result = cas9_view(df, editor)
    else:
        df = BE
        if omim:
            # Convert omim to float for comparison, handles both "136880" and "136880.0"
            try:
                omim_val = float(omim)
                df = df[df["OMIM_ID"] == omim_val]
            except:
                df = df[df["OMIM_ID"].astype(str)==omim]
        if mutation:
            # case-insensitive substring match on Name
            df = df[df["Name"].str.contains(mutation, case=False, na=False, regex=False)]
        result = be_view(df, editor)

    # sorting/grouping
    if mode == "be":
        if sort_by and sort_by in result.columns:
            # Sort by selected column descending only
            result = result.sort_values(by=sort_by, ascending=False)
        elif "Prism-CRISPRScore" in result.columns:
            # Default sort by Prism-CRISPRScore descending
            result = result.sort_values(by="Prism-CRISPRScore", ascending=False)
    else:
        # default sort: group by Name then Prism-CRISPRScore desc
        # if custom sort selected: sort by that column descending only
        if sort_by and sort_by in result.columns:
            result = result.sort_values(by=sort_by, ascending=False)
        elif "Name" in result.columns and "Prism-CRISPRScore" in result.columns:
            result = result.sort_values(by=["Name","Prism-CRISPRScore"], ascending=[True,False])

    # ecords = dataframe_to_ordered_records(result.head(200))
    records = dataframe_to_ordered_records(result)
    json_str = json.dumps(records, separators=(',', ': '))
    return Response(json_str, mimetype='application/json')

@app.route("/api/search_original", methods=["POST"])
def search_original():

    try:

        load_all_data()
        
        data = request.json
        omim = data.get("omim_id","").strip()
        mutation = data.get("mutation","").strip()
        mode = data.get("mode")
        editor = data.get("editor")
        name = data.get("name")
        sgRNA = data.get("sgRNA")
        sort_by = data.get("sort_by")  # 添加sort_by参数
        
        if mode != "be" or editor is None:
            return jsonify([])

        global GR_ABE, GR_CBE, DEEPBE_ABE_SCORES, DEEPBE_CBE_SCORES, BEDEE_ABE_SCORES, BEDEE_CBE_SCORES
        
        abe = GR_ABE.copy()
        cbe = GR_CBE.copy()
        df = pd.concat([abe,cbe], ignore_index=True)
        
        df["target_base"] = df["sgrna"].astype(str).str[:20].str.upper()
        df["outcome_base"] = df["outcome"].astype(str).str[:20].str.upper()
        
        if omim:
            # Convert omim to float for comparison, handles both "136880" and "136880.0"
            try:
                omim_val = float(omim)
                df = df[df["OMIM_ID"] == omim_val]
            except:
                df = df[df["OMIM_ID"].astype(str)==omim]
        if mutation:
            # case-insensitive substring match on Name
            df = df[df["Name"].str.contains(mutation, case=False, na=False, regex=False)]
        if name:
            # exact match on Name
            df = df[df["Name"] == name]
        if sgRNA:
            # exact match on target_base (sgRNA)
            df = df[df["target_base"] == sgRNA]
        
        df = df.merge(DEEPBE_ABE_SCORES, on=["Name", "target_base", "outcome_base"], how="left")
        df = df.merge(DEEPBE_CBE_SCORES, on=["Name", "target_base", "outcome_base"], how="left")
        df = df.merge(BEDEE_ABE_SCORES, on=["Name", "target_base", "outcome_base"], how="left")
        df = df.merge(BEDEE_CBE_SCORES, on=["Name", "target_base", "outcome_base"], how="left")
        
        result = be_view_original(df, editor)

        # sorting
        if sort_by and sort_by in result.columns:
            # Sort by selected column descending only
            result = result.sort_values(by=sort_by, ascending=False)
        elif "Prism-CRISPRScore" in result.columns:
            # Default sort by Prism-CRISPRScore descending
            result = result.sort_values(by="Prism-CRISPRScore", ascending=False)

        # records = dataframe_to_ordered_records(result.head(200))
        records = dataframe_to_ordered_records(result)
        json_str = json.dumps(records, separators=(',', ': '))
        return Response(json_str, mimetype='application/json')
    except Exception as e:
        print(f"Error in search_original: {e}")
        return jsonify({"error": str(e)}), 500


def cas9_view(df, editor):
    mapping = {
        "WT-SpCas9": ("Prism-CRISPR_WT","DeepMEns_WT","TransCrispr_WT"),
        "eSpCas9(1.1)": ("Prism-CRISPR_ESP","DeepMEns_ESP","TransCrispr_ESP"),
        "SpCas9-HF1": ("Prism-CRISPR_HF","DeepMEns_HF","TransCrispr_HF")
    }
    gr, deep, trans = mapping[editor]
    result = pd.DataFrame({
        "Name": df["Name"],
        "OMIM_ID": pd.to_numeric(df["OMIM_ID"], errors='coerce').astype('Int64'),
        "sgRNA": df["sgRNA_base"],
        "PAM": df["pam"],
        "Prism-CRISPRScore": df[gr],
        "DeepMEnsScore": df[deep],
        "TransCrisprScore": df[trans],
        "Editor": editor
    })
    return result[["Name", "OMIM_ID", "sgRNA", "PAM", "Prism-CRISPRScore", "DeepMEnsScore", "TransCrisprScore", "Editor"]]

def be_view(df, editor):
    mapping = {
        "ABE7.10": ("Prism-CRISPR_ABE7","DeepBE_ABE7","BEDeepon_ABE7"),
        "ABEmax": ("Prism-CRISPR_ABEmax","DeepBE_ABEmax","BEDeepon_ABEmax"),
        "ABE8e": ("Prism-CRISPR_ABE8e","DeepBE_ABE8e","BEDeepon_ABE8e"),
        "BE4": ("Prism-CRISPR_BE4","DeepBE_BE4","BEDeepon_BE4"),
        "CBE4max": ("Prism-CRISPR_CBE4max","DeepBE_CBE4max","BEDeepon_CBE4max"),
        "Target-AID": ("Prism-CRISPR_AID","DeepBE_AID","BEDeepon_AID")
    }
    gr, deep, deepon = mapping[editor]
    
    df = df.copy()
    df["sgRNA_trimmed"] = df["target_base"].str[:20].str.upper()
    df["outcome_trimmed"] = df["outcome"].astype(str).str[:20].str.upper()
    
    df = df[df["sgRNA_trimmed"] == df["outcome_trimmed"]]
    
    result = pd.DataFrame({
        "Name": df["Name"],
        "OMIM_ID": pd.to_numeric(df["OMIM_ID"], errors='coerce').astype('Int64'),
        "sgRNA": df["target_base"],
        "PAM": df["pam"],
        "Outcome": df["outcome"].astype(str).str[:20],
        "Prism-CRISPRScore": 1 - df[gr],
        "DeepBEScore": 1 - df[deep],
        "BEDeeponScore": 1 - df[deepon],
        "Editor": editor
    })
    
    # Drop rows with NaN scores
    result = result.dropna(subset=["Prism-CRISPRScore", "DeepBEScore", "BEDeeponScore"])
    
    result = result.drop_duplicates(subset=["Name", "sgRNA"], keep="first")
    
    return result[["Name", "OMIM_ID", "sgRNA", "PAM", "Outcome", "Prism-CRISPRScore", "DeepBEScore", "BEDeeponScore", "Editor"]]

def be_view_original(df, editor):
    mapping = {
        "ABE7.10": ("Prism-CRISPR_ABE7","DeepBE_ABE7","BEDeepon_ABE7"),
        "ABEmax": ("Prism-CRISPR_ABEmax","DeepBE_ABEmax","BEDeepon_ABEmax"),
        "ABE8e": ("Prism-CRISPR_ABE8e","DeepBE_ABE8e","BEDeepon_ABE8e"),
        "BE4": ("Prism-CRISPR_BE4","DeepBE_BE4","BEDeepon_BE4"),
        "CBE4max": ("Prism-CRISPR_CBE4max","DeepBE_CBE4max","BEDeepon_CBE4max"),
        "Target-AID": ("Prism-CRISPR_AID","DeepBE_AID","BEDeepon_AID")
    }
    gr, deep, deepon = mapping[editor]
    
    result = pd.DataFrame({
        "Name": df["Name"],
        "OMIM_ID": pd.to_numeric(df["OMIM_ID"], errors='coerce').astype('Int64'),
        "sgRNA": df["target_base"],
        "PAM": df["pam"],
        "Outcome": df["outcome"].astype(str).str[:20],
        "Prism-CRISPRScore": df[gr], 
        "DeepBEScore": df[deep],    
        "BEDeeponScore": df[deepon], 
        "Editor": editor
    })
    
    # Drop rows with NaN scores
    result = result.dropna(subset=["Prism-CRISPRScore", "DeepBEScore", "BEDeeponScore"])
    
    return result[["Name", "OMIM_ID", "sgRNA", "PAM", "Outcome", "Prism-CRISPRScore", "DeepBEScore", "BEDeeponScore", "Editor"]]

app = app

if __name__=="__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)