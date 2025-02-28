# **Project name: Type 1 Diabetes Software**
***
The software is a web-based application that allows users to search for data related to single nucleotide polymorphisms (SNPs) on chromosome 6 that are associated with increased risk of Type 1 Diabetes (T1D).
***
## **Table of contents**
***
1. General Info
2. Software Technologies
3. Data Sources
4. Installation
5. Methodology
6. FAQ's

### **General Info**
***
The web application is based off of Python 3 making extensive use of the Flask package to generate HTML web-pages from the Python script, this communicates with the SQLite database to retrieve and store user data and processes user requests, generating an appropriate HTML response using the templates. When you use the web application you will be prompted to input one of three things: an rsID, a gene name or a location range on chromosome 6. After inputting one of these values the following results will be shown. If a singular SNP is returned from the search the user will see the following information: SNP name, genomic position, p-value from the association test, mapped gene, variant frequency in three different human populations of interest, one measure of functional impact and/or clinical relevance for each variant, and one functional or gene ontology term associated with each mapped gene. If multiple SNP’s are returned from your search in addition to this information, a Manhattan Plot, Linkage Disequilibrium plot and a text file of LD values will be presented.

### **Libraries for Installation**
***

* SQLite3 
* numpy
* flask
* pandas
* matplotlib
* wtforms
* requests
* ld_plot
* pathlib
***

### **Methodology**
* Install all the latest version of the required libraries listed above from python using the following script:
   pip install -r requirements.txt

* When making the database, ensure you have both CSV file in the same directory before running database_code.py 

* In order to successfully run the web application, ensure the database file and Web_app_final.py are in the same directory and that the directory contains a subdirectory named templates, containing all of the HTML files

* Using Jupyter, open a new terminal, change directory to the directory containing the Web_app_final.py file, and enter the following script to run the software:
   Python Web_app_final.py

* An input query of rsID returns all the available information related to that in the entry search for chromosome 6.

* An input of a gene name will return information for all SNP's associated with that gene. 

* An input of location range for chromosome 6 returns information relating to all SNPs located within that range, a Manhattan plot, and an option to choose two or more SNPs for which to calculate downloadable linkage disequilibrium (LD) scores and produce an LD plot. 
***

## **FAQ's**
* **I ran pip install -r requirements.txt but got an error opening the mentioned file in the terminal**

_If the command fails, a common cause of this is your terminal is running the command in the same directory as requirments.txt_

* **I give an rsid in the search entry and theres no result displayed**

_Avoid blank spaces after inputting an rsID, also make sure to use an rsID or a gene name that is relevant to Type 1 Diabetes_

* **I have searched chromosome positions that are related to T1D but nothing comes up**

_This web app is specific to chromosome 6 only, for now please keep your searches limited to this chromosome_
