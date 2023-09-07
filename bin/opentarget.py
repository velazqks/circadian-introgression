#!/usr/bin/env python
# opentarget.py
"""""""""""""""""""""""""""""""""
# Author: Sarah Fong
#
# Description: Import relevant libraries for HTTP request and JSON formatting
#
"""""""""""""""""""""""""""""""""


import requests
import json
import pandas as pd

def query_OpenTargetGenetics_genesForVariant(variant_id):
    
    # Build query string
    query_string = """
    query genesForVariant ($myVariantId: String!){
      genesForVariant(variantId:$myVariantId) {
        gene {
            id
            symbol
            description
            }
        variant
        overallScore
        qtls {
            typeId
            sourceId
            aggregatedScore
            tissues {
                tissue {
                    id
                    name
                    }
                }
                
            }
        distances {
        tissues { distance
        }
        }
        
    }
}
    """ 

    # Set variables object of arguments to be passed to endpoint
    variables = { "myVariantId": variant_id}

    # Set base URL of Genetics Portal GraphQL API endpoint
    base_url = "https://api.genetics.opentargets.org/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})
    #print(r.status_code)

    # Transform API response into JSON 
    api_response_as_json = json.loads(r.text)
    return api_response_as_json
    """
    # Length of the results 
    #print(api_response_as_json)
    #res_len = len(api_response_as_json["data"]["genesForVariant"]["gene".values()])

    # Return results or flag of no results
    #if res_len >0:
        
        # Print first element of JSON response data
        #print(api_response_as_json["data"]["pheWAS"]["associations"][0])
        
        # Return results
        return api_response_as_json
    
    else:
        # Flag no results
        
        #print("no results")
        
        # Return None
        return None
    """

    
    
def gene_outputs_todict(query_results):
    
    """
    parse query_OpenTargetGenetics_genesForVariant results for one single variant. 
    
    1 - make dictionary to collect results
    2 - per key (field/measurement), add value to list of values
    3 - if key is the study object, the study object is a dictionary and we have to parse this separately
        3a - if study information is None, rename j,k variables as key = None, value = empty_dict
    4 - if key is not already in dictionary, create key w empty list, append value to list, then reassign dictionary key w/ value list
    5 - if key is in dictionary, append value to list, then reassign dictionary key w/ value list
    
    """
    
    #1
    results_dict = {}
    
    #2
    for i in query_results["data"]["genesForVariant"]["gene"]:    
        for key, val in i.items():
            
            #3
            if type(key) == dict:
                
                #3a
                if val == None:
                    continue
                    
                for j,k in val.items():
                    
                    #4
                    if j not in results_dict.keys():
                        results_dict[j] = []
                        results_list = results_dict[j]
                        results_list.append(k)
                        results_dict[j] = results_list
                    #5
                    else:
                        results_list = results_dict[j]
                        results_list.append(k)
                        results_dict[j] = results_list

            else:
                
                #4
                if key not in results_dict.keys():
                    results_dict[key] = []
                    results_list = results_dict[key]
                    results_list.append(val)
                    results_dict[key] = results_list
                #5
                else:
                    results_list = results_dict[key]
                    results_list.append(val)
                    results_dict[key] = results_list
    return results_dict


###
# make lists of genes, eQTL genes, and genes within 500kb
###

def get_genes_per_var(ID):
    
    info = query_OpenTargetGenetics_genesForVariant(ID)

    ## make lists of genes, eQTL genes, and genes within 500kb

    genes, eqtl_genes, megabase = [],[],[]

    for v in info.values():
        for i in v.values():
            for j in i:

                gene = j["gene"]['symbol']
                eqtl = j["qtls"]
                dis = j["distances"]
                
                if len(eqtl)>0:
                    eqtl_genes.append(gene)

                if len(dis)>0:
                    distance = dis[0]["tissues"][0]["distance"]
                
                    if distance< 5e5:
                        megabase.append(gene)

                genes.append(gene)
                
    return genes, eqtl_genes, megabase
###
# PHEWAS
###
# Set study_id and variant_id variables
# variant_id = "11_72752390_G_A"

def query_OpenTargetGenetics_pheWAS(variant_id):
    
    # Build query string
    query_string = """
        query pheWAS ($myVariantId: String!) {
          pheWAS(variantId:$myVariantId) {
            associations {
              studyId
              eaf
              beta
              se
              nTotal
              nCases
              oddsRatio
              pval
              study {
                traitCategory
                traitEfos 
                traitReported
                numAssocLoci
              }
            }
          }
        }

    """ 
    # Efo = "experimental factor ontology"

    # Set variables object of arguments to be passed to endpoint
    variables = { "myVariantId": variant_id}

    # Set base URL of Genetics Portal GraphQL API endpoint
    base_url = "https://api.genetics.opentargets.org/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})
    #print(r.status_code)

    # Transform API response into JSON 
    api_response_as_json = json.loads(r.text)

    # Length of the results 
    res_len = len(api_response_as_json["data"]["pheWAS"]["associations"])

    # Return results or flag of no results
    if res_len >0:
        
        # Print first element of JSON response data
        #print(api_response_as_json["data"]["pheWAS"]["associations"][0])
        
        # Return results
        return api_response_as_json
    
    else:
        # Flag no results
        
        #print("no results")
        
        # Return None
        return None

###
# turn open target PHEWAS outputs into a dictionary
###

def pheWAS_outputs_todict(query_results):
    
    """
    parse the open target results for one single variant. 
    
    1 - make dictionary to collect results
    2 - per key (field/measurement), add value to list of values
    3 - if key is the study object, the study object is a dictionary and we have to parse this separately
        3a - if study information is None, rename j,k variables as key = None, value = empty_dict
    4 - if key is not already in dictionary, create key w empty list, append value to list, then reassign dictionary key w/ value list
    5 - if key is in dictionary, append value to list, then reassign dictionary key w/ value list
    
    """
    
    #1
    results_dict = {}
    
    #2
    for i in query_results["data"]["pheWAS"]["associations"]:    
        for key, val in i.items():
            
            #3
            if key == "study":
                
                #3a
                if val == None:
                    continue
                    
                for j,k in val.items():
                    
                    #4
                    if j not in results_dict.keys():
                        results_dict[j] = []
                        results_list = results_dict[j]
                        results_list.append(k)
                        results_dict[j] = results_list
                    #5
                    else:
                        results_list = results_dict[j]
                        results_list.append(k)
                        results_dict[j] = results_list

            else:
                
                #4
                if key not in results_dict.keys():
                    results_dict[key] = []
                    results_list = results_dict[key]
                    results_list.append(val)
                    results_dict[key] = results_list
                #5
                else:
                    results_list = results_dict[key]
                    results_list.append(val)
                    results_dict[key] = results_list
    return results_dict

###
# turn PheWAS dictionary into a dataframe 
###

def dict_to_df(results_dict, var_id):

    """
    1 - create object newdf
    2 - parse through results_dict
    
    3 - if this is the first val, make a dataframe of the new df
    4 - in the key is "nCases", this key contains duplicate elements, so I selected every other element
    5 - append key, values as new column and rows
    """
    newdf = None
    n_val = 0
    
    for key, val in results_dict.items():
        val_len = len(val)
        
        if newdf is None:
            newdf = pd.DataFrame({key:val})
            n_val = val_len
    
        if val_len == n_val:
            newdf[key] = val
        
        elif key == "nCases":  # nCases has duplicate values, so pick every other 
            
            newdf[key] = val[::1]

    newdf["var_id"] = var_id
                
    return newdf

###
# get PheWAS from the open genetics pipeline
###

def run_otg_phewas(var_list):
     
    """
    return pandas dataframe of variants and pheWAS overlaps from UKBB via open targets genetics. 
    
    input
        variant list (list) where each variant is a string like this: 1_154453788_C_T or CHR_POS_REF_ALT
    
    output
        pandas dataframe (dataframe) of pheWAS variants and traits that overlap variant. 
        
    method
    
        1 - dictionaries for collecting open target results + atac-starr feature results
        2 - per variant
        2a - query open target genetics for phewas results linked to that variant. 
        3 - if pheWAS results, turn results into a dictionary for that variant
        4 - if variant dictionary, turn dictionary into a dataframe
        5 - save each variant open target and huacc/atac-starr results to a dictionary 
        6 - if there are results, concatenate variant results and return dataframe. 
        
    """
    
    #1
    res = {}  # collect open target results
    
    #2
    for var_id in var_list:
        
        
        #2a
        query_results = query_OpenTargetGenetics_pheWAS(var_id)  # query opentargetgenetics
        
        #3
        if query_results != None and var_id != None:
            if var_id not in res.keys():
                
                #3
                results_dict = pheWAS_outputs_todict(query_results)  # turn results into dictionary
                
                #4
                var_df = dict_to_df(results_dict, var_id)  # turn pheWAS outputs dictionary into a dataframe
                
                #5
                res[var_id] = var_df  # add results to dictionary
                
            else:
                continue
                #print("tested this var_id already")
        else: 
            continue
            #print("no pheWAS")
    # 6
    if len(res.keys())>0:
        
        otg_phewas_results = pd.concat(res.values())  # open target phewas results 
        

    else:
        otg_phewas_results = None
    
    return otg_phewas_results
