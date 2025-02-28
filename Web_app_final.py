# Importnat necessary packges
from flask import Flask, render_template, url_for, redirect, request, make_response
import pandas as pd
import sqlite3
import io
import base64
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.pyplot as plt
import urllib.request

from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import InputRequired

import requests
import re
import numpy as np
from ld_plot.ld_plot import ld_plot
import matplotlib.pyplot as plt
import itertools

# Create the Flask instance and a secret key
app = Flask(__name__)
app.config['SECRET_KEY'] = 'Firefox'

# Create a connection to the database 
con = sqlite3.connect("app_data_comp.db", check_same_thread=False)
cur = con.cursor()
con.row_factory = sqlite3.Row

# Create classes for each user submission box
class rsIDForm(FlaskForm):
    rs_ID = StringField('Enter a valid rsID value:', validators=[InputRequired()])
    submit_rsID = SubmitField('Search')
    
class GeneForm(FlaskForm):
    gene_name = StringField('Enter a valid gene name:', validators=[InputRequired()])
    submit_gene = SubmitField('Search')

class LocationForm(FlaskForm):
    chromosome = StringField('Chromosome:', validators=[InputRequired()])
    startbp = StringField('Start basepair (bp):', validators=[InputRequired()])
    endbp = StringField('End basepair (bp):', validators=[InputRequired()])
    submit_location = SubmitField('Search')

class QueryForm(FlaskForm):
    rsvalue = StringField('Insert space seperated rsID values to calculate linkage disequilibrium:', validators=[InputRequired()]) 
    submit_LD = SubmitField('Calculate Linkage Disequilibrium')

# Define the action for the top level route
@app.route('/')
def index():
    return render_template('base_index.html')

# Define action for SNP route to retrieve data based on rsID input
@app.route('/SNP', methods=['GET', 'POST'])
def SNP():
    # Assign response of rsIDForm to variable 
    rform = rsIDForm()
    
    # Validation of user input and data extraction
    if rform.validate_on_submit():
        
        # Get rs_ID value from form data
        rs_ID = rform.rs_ID.data.lower()
        
        # Query the database
        data = pd.read_sql_query('''SELECT C6_data.rsID, location, risk_allele, pvalue, variant_functional_impact, \
               variant_clinical_relevance, FIN_freq, JPT_freq, NGA_freq, mapped_gene_1, mapped_gene_2, \
               biological_process_1, cellular_component_1, molecular_function_1, biological_process_2, cellular_component_2, \
               molecular_function_2 FROM C6_data, GO_terms WHERE C6_data.rsID = GO_terms.rsID''', con, index_col='rsID')
        
        # Convert database query output to dataframe
        df = pd.DataFrame(data, columns=["location", "rsID", "risk_allele", "pvalue", "variant_functional_impact", \
                "variant_clinical_relevance", "FIN_freq", "JPT_freq", "NGA_freq", "mapped_gene_1", "mapped_gene_2", \
                "biological_process_1", "cellular_component_1", "molecular_function_1", "biological_process_2", "cellular_component_2", \
                "molecular_function_2"])
       
        try:
            # Locate user search query from dataframe 
            row = df.loc[rs_ID]
            
            # Rendering of data extraction in a new webpage
            return render_template('rsID.html', ID=rs_ID, Location=row["location"], Allele=row["risk_allele"], pval=row["pvalue"], \
                   VFI=row["variant_functional_impact"], VCR=row["variant_clinical_relevance"], FINf=row["FIN_freq"], \
                   JPTf=row["JPT_freq"], NGAf=row["NGA_freq"], Gene1=row["mapped_gene_1"], \
                   BP1=row["biological_process_1"], CC1=row["cellular_component_1"], MF1=row["molecular_function_1"], \
                   Gene2=row["mapped_gene_2"], BP2=row["biological_process_2"], \
                   CC2=row["cellular_component_2"], MF2=row["molecular_function_2"])    
        except:
            # If rsID is not found the below error message is received
            return "Information unavailable for rsID %s." % rs_ID
    
    # Render web page for rsID submission box
    return render_template('r_index.html', rform=rform)

# Define action for Gene route to retrieve data based on Gene input
@app.route('/Gene', methods=['GET', 'POST'])
def Gene():
    # Assign response of GeneForm to variable 
    gform = GeneForm()
    
    # Validation of user input and data extraction
    if gform.validate_on_submit():
        
        # Get gene_ value from form data
        gene_name = gform.gene_name.data.upper() 
        
        # Query the database
        data = pd.read_sql_query("""SELECT C6_data.rsID, chromosome, location, risk_allele, pvalue, \
                variant_functional_impact, variant_clinical_relevance, FIN_freq, JPT_freq, NGA_freq, \
                mapped_gene_1, mapped_gene_2, biological_process_1, cellular_component_1, molecular_function_1, \
                biological_process_2, cellular_component_2, molecular_function_2 \
                FROM C6_data, GO_terms WHERE C6_data.rsID=GO_terms.rsID""", con)
        
        # Convert database query output to dataframe
        table = pd.DataFrame(data, columns=["chromosome", "location", "rsID", "risk_allele", "pvalue", \
                "variant_functional_impact", "variant_clinical_relevance", "FIN_freq", "JPT_freq", \
                "NGA_freq", "mapped_gene_1", "mapped_gene_2", "biological_process_1", "cellular_component_1", \
                "molecular_function_1", "biological_process_2", "cellular_component_2", \
                "molecular_function_2"])
        try:
            # Locate user search query from dataframe 
            row = table[table.eq(gene_name).any(1)]
            
            # Rendering of data extraction in a new webpage
            return render_template('g_table.html', Gene = gene_name, tables=[row.to_html(classes='data')], titles=row.columns.values)
 
        except:
            # If gene is not found the below error message is received
            return "Information unavailable for gene %s." % gene_name
        
    # Render web page for gene name submission box
    return render_template('g_index.html', gform=gform)    
    
# Define action for Location route to retrieve data based on Location input
@app.route('/Location', methods=['GET', 'POST'])
def Location():
    # Assign response of LocationForm to variable 
    lform = LocationForm()
    
    # Validation of user input and data extraction
    if lform.validate_on_submit():
        chromosome = lform.chromosome.data
        startbp = lform.startbp.data
        endbp = lform.endbp.data
        
        # Query the database
        data = pd.read_sql_query("""SELECT C6_data.rsID, chromosome, location, risk_allele, pvalue, \
            variant_functional_impact, variant_clinical_relevance, FIN_freq, JPT_freq, NGA_freq, \
            mapped_gene_1, mapped_gene_2, cellular_component_1, molecular_function_1, \
            biological_process_1, cellular_component_2, molecular_function_2, biological_process_2 \
            FROM C6_data, GO_terms WHERE C6_data.rsID = GO_terms.rsID AND chromosome=%s AND location BETWEEN %s \
            AND %s """ % (chromosome, startbp, endbp), con, index_col='location')
        
        # Convert database query output to dataframe
        table = pd.DataFrame(data, columns=["chromosome", "rsID", "risk_allele", "pvalue", \
                "variant_functional_impact", "variant_clinical_relevance", "FIN_freq", "JPT_freq", \
                "NGA_freq", "mapped_gene_1", "mapped_gene_2", "cellular_component_1", \
                "molecular_function_1", "biological_process_1", "cellular_component_2", \
                "molecular_function_2", "biological_process_2"])
        
        # Second database query to extract data for the Manhattan plot
        queryM = f"SELECT * FROM C6_data WHERE chromosome = '{chromosome}' AND location BETWEEN {startbp} AND {endbp}"
        df = pd.read_sql(queryM, con)
        df['pvalue'] = df['pvalue'].apply(lambda x: float(x))
    
        # Group the data by chromosome and find the -log10(p-value) for each SNP
        grouped = df.groupby('chromosome')
        x = []
        y = []
        labels = []
        for name, group in grouped:
            labels += list(group['rsID'])
            x += list(group['location'])
            y += [-1 * np.log10(float(pvalue)) if float(pvalue) > 0 else np.nan for pvalue in group['pvalue']]
    
        # Generate Manhattan plot if multiple SNPs are present
        if len(x) > 1:
            fig, ax = plt.subplots()
            ax.scatter(x, y, s=3)
            ax.set_xlabel('Base pair position')
            ax.set_ylabel('-log10(p-value)')
            ax.set_title('Manhattan plot')
            img=io.BytesIO()
            plt.savefig(img, format='png')
            img.seek(0)
            plot_data = urllib.parse.quote(base64.b64encode(img.getvalue()).decode('utf-8'))

        try:
            # Rendering of data extraction in a new webpage including Manhattan plot
            return render_template('l_table.html',  tables=[table.to_html(classes='data')], \
                                   titles=table.columns.values, plot_url=plot_data)

        except:
            # If no rsIDs found in location range the below error message is received
            return "No Type 1 Diabetes associated variants have been found in genomic coordinates %s:%s-%s." % (chromosome, startbp, endbp)
    
    # Render web page for genomic coordinates submission boxes
    return render_template('l_index.html', lform=lform)

# Define action for R2_index route 
@app.route('/R2_index', methods=['GET','POST'])
def R2_index():
    
    # Assign response of QueryForm to variable 
    LDform = QueryForm()
    
    # Validation of user input, data exrtaction, and redirection to data output web page
    if LDform.validate_on_submit():
        rsvalue = LDform.rsvalue.data
        return redirect(url_for('R2_value', rsvalue=rsvalue))
    
    # Render web page for rsID submission boxes
    return render_template('LD_index.html', LDform=LDform)

# App route for function to extract R2 values into seperate lists and dataframes for each country
@app.route('/R2_value/<rsvalue>')    
def R2_value(rsvalue):
    
    # Split the input to form a list
    rsIDs = rsvalue.split()
        
    # Seperate rsIDs into all possible combinations    
    rsID_combinations = [(rsIDs[i], rsIDs[j]) for i in range(len(rsIDs)) for j in range(i+1, len(rsIDs))]
    
    # Extract relevant R2 values into lists
    R2 = []
    for x in rsID_combinations:
            
        rsID1 = x[0]
        rsID2 = x[1]

        # Create variable with API request information
        LDlink_request = f'https://ldlink.nci.nih.gov/LDlinkRest/ldpop?var1={rsID1}&var2={rsID2}&pop=FIN%2BJPT%2BYRI&r2_d=r2&genome_build=grch38&token=4c5ec6b51d7e'
            
        # API request to LDlink website for LDpop
        r = requests.get(LDlink_request)

        # Read output into string and format correctly
        txt = r.text
        txt_split = re.split(r"\t|\n", txt)
        txt_split.pop()

        # Split list into four lists for each row of information
        split_list = []
        for k in range(0, len(txt_split), 9):
                split_list.append(txt_split[k:k + 9])

        # Remove list of headers
        split_list.pop(0)

        # Extract R2 values and multiply by 100
        R2_temp = []
        for p in split_list:
            R2_temp.append(p[5])
        for d in range(len(R2_temp)):
            if R2_temp[d] == "NA":
                       R2_temp[d] = 0
            R2_temp[d] = float(R2_temp[d]) * 100
        R2.append(R2_temp)    

    # Split output R2 list into 3 seperate lists, one for for each country
    NGA_R2 = []
    JPT_R2 = []
    FIN_R2 = []
    
    for z in R2:
        NGA_R2.append(z[0])
        JPT_R2.append(z[1])
        FIN_R2.append(z[2])    
    
    # Make symmetric matrix for each country

    sm_NGA = np.zeros([len(rsIDs),len(rsIDs)], dtype=np.double)
    x1,y1 = np.triu_indices(len(rsIDs), k=1)
    sm_NGA[x1,y1] = NGA_R2
    sm_NGA[y1,x1] = NGA_R2
    sm_NGA[ np.diag_indices(len(rsIDs)) ] = 100
       
    sm_JPT = np.zeros([len(rsIDs),len(rsIDs)], dtype=np.double)
    x2,y2 = np.triu_indices(len(rsIDs), k=1)
    sm_JPT[x2,y2] = JPT_R2
    sm_JPT[y2,x2] = JPT_R2
    sm_JPT[ np.diag_indices(len(rsIDs)) ] = 100
    
    sm_FIN = np.zeros([len(rsIDs),len(rsIDs)], dtype=np.double)
    x3,y3 = np.triu_indices(len(rsIDs), k=1)
    sm_FIN[x3,y3] = FIN_R2
    sm_FIN[y3,x3] = FIN_R2
    sm_FIN[ np.diag_indices(len(rsIDs)) ] = 100
    
     # Split NGA matrix into lists for dataframe input
    NGA_list = []
    for m_N in sm_NGA:
        NGA_list.append(list(m_N))
        
    # Make dataframe of R2 values for NGA
    rsID_index = rsIDs  
    
    #global NGA_df
    NGA_df = pd.DataFrame(index = rsID_index, columns = rsIDs)
    
    col = 0   
    row = 0   
    
    for a_N in NGA_list:
        for b_N in a_N:
            NGA_df.at[rsID_index[row], rsIDs[col]] = b_N
            col +=1
            if col == len(rsIDs):
                col = 0
        row += 1
    
    # Split JPT matrix into lists for dataframe input
    JPT_list = []
    for m_J in sm_JPT:
        JPT_list.append(list(m_J))
        
    # Make dataframe of R2 values for JPT 
    global JPT_df
    JPT_df = pd.DataFrame(index = rsID_index, columns = rsIDs)
    
    col = 0   
    row = 0   
    
    for a_J in JPT_list:
        for b_J in a_J:
            JPT_df.at[rsID_index[row], rsIDs[col]] = b_J
            col +=1
            if col == len(rsIDs):
                col = 0
        row += 1
    
    # Split FIN matrix into lists for dataframe input
    FIN_list = []
    for m_F in sm_FIN:
        FIN_list.append(list(m_F))
        
    # Make dataframe of R2 values for FIN
    global FIN_df
    FIN_df = pd.DataFrame(index = rsID_index, columns = rsIDs)
    
    row = 0
    col = 0
    
    for a_F in FIN_list:
        for b_F in a_F:
            FIN_df.at[rsID_index[row], rsIDs[col]] = b_F
            col +=1
            if col == len(rsIDs):
                col = 0
        row += 1
    
    # Make plot for NGA
    NGA_plot = ld_plot(ld=sm_NGA, labels=rsIDs)
    plt.title("Linkage disequilibrium plot for Nigeria (NGA)", loc = "left")
    plt.xticks(fontsize=7.5)
    img_NGA=io.BytesIO()
    plt.savefig(img_NGA, format='png')
    img_NGA.seek(0)
    plot_data_NGA = urllib.parse.quote(base64.b64encode(img_NGA.getvalue()).decode('utf-8'))
    
    # Make plot for JPT
    JPT_plot = ld_plot(ld=sm_JPT, labels=rsIDs)
    plt.title("Linkage disequilibrium plot for Japan (JPT)", loc = "left")
    plt.xticks(fontsize=7.5)
    img_JPT=io.BytesIO()
    plt.savefig(img_JPT, format='png')
    img_JPT.seek(0)
    plot_data_JPT = urllib.parse.quote(base64.b64encode(img_JPT.getvalue()).decode('utf-8'))
    
    # Make plot for FIN
    FIN_plot = ld_plot(ld=sm_FIN, labels=rsIDs)
    plt.title("Linkage disequilibrium plot for Finland (FIN)", loc = "left")
    plt.xticks(fontsize=7.5)
    img_FIN=io.BytesIO()
    plt.savefig(img_FIN, format='png')
    img_FIN.seek(0)
    plot_data_FIN = urllib.parse.quote(base64.b64encode(img_FIN.getvalue()).decode('utf-8'))
    
    # Return the tables and plots
    return render_template('LD_graphs.html', \
           NGA_table=[NGA_df.to_html(classes='NGA_df').replace('\n', '')], NGA_titles=NGA_df.columns.values, \
           JPT_table=[JPT_df.to_html(classes='JPT_df').replace('\n', '')], JPT_titles=JPT_df.columns.values, \
           FIN_table=[FIN_df.to_html(classes='FIN_df').replace('\n', '')], FIN_titles=FIN_df.columns.values, plot_url_NGA=plot_data_NGA, plot_url_JPT=plot_data_JPT, plot_url_FIN=plot_data_FIN)

@app.route('/download_NGA')
def download_NGA():
    resp = make_response(NGA_df.to_csv())
    resp.headers["Content-Disposition"] = "attachment; filename=Nigeria_LD_data.csv"
    resp.headers["Content-Type"] = "text/csv"
    return resp

@app.route('/download_JPT')
def download_JPT():
    resp = make_response(JPT_df.to_csv())
    resp.headers["Content-Disposition"] = "attachment; filename=Japan_LD_data.csv"
    resp.headers["Content-Type"] = "text/csv"
    return resp

@app.route('/download_FIN')
def download_FIN():
    resp = make_response(FIN_df.to_csv())
    resp.headers["Content-Disposition"] = "attachment; filename=Finland_LD_data.csv"
    resp.headers["Content-Type"] = "text/csv"
    return resp

if __name__ == '__main__':
    app.run(debug=True)