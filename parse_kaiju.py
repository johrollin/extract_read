#!/usr/bin/env python3

import argparse
import pandas as pd
from ete3 import NCBITaxa
import json

def get_kaiju_data(result_file_path):
    """
    Open tab kaiju result file into dataframe
    store unclassified reads id in df_U_data
    store classified reads id in df_U_data (with tax_id and score ...)
    """
    pd.set_option('mode.chained_assignment', None) # !!!! careful this modification may hide error if df_all_data is use in the futur 
    data_file = open(result_file_path, "r")
    try:
        df_all_data = pd.read_csv(data_file, sep="\t",dtype="str", \
        names=["Classified", "Read_name", "ncbi_id", "score", "ncbi_id_best_hit", \
        "accession_num", "matching_seq"])
        df_U_data = df_all_data[~df_all_data.Classified.str.contains("C")] # store unclassified reads
        df_U_data.drop(["ncbi_id", "score", "ncbi_id_best_hit","accession_num", "matching_seq"]\
        ,axis=1, inplace=True ) 
        df_C_data = df_all_data[~df_all_data.Classified.str.contains("U")] # remove unclassified reads
        #df_C_data.drop(["ncbi_id_best_hit", "matching_seq"],axis=1, inplace=True ) # remove unwanted column to save RAM
    except pandas.errors.ParserError: ###  pas sur que l'erreur se represente 
        df_all_data = pd.read_csv(data_file, sep="\t", dtype="str", \
        names=["Classified", "Read_name", "ncbi_id"])
        df_U_data = df_all_data[~df_all_data.Classified.str.contains("C")] # store unclassified reads
        df_C_data = df_all_data[~df_all_data.Classified.str.contains("U")] # remove unclassified reads
    if df_C_data["accession_num"].isnull().values.any():
        prot_analysis_possible=False
    else:
        prot_analysis_possible=True
    reads_tot = len(df_U_data) + len(df_C_data)
    msg_U = "NB of Unclassified reads: " + str(len(df_U_data)) \
        + " " + str(round((len(df_U_data)/reads_tot*100), 2)) + "%"
    msg_C = "NB of Classified reads: " + str(len(df_C_data)) \
        + " " + str(round((len(df_C_data)/reads_tot)*100, 2)) + "%"
    print(msg_U)
    print(msg_C)
    data_file.close()

    return df_C_data, df_U_data, reads_tot, prot_analysis_possible

def explore_protein_data(df_protein, filter_term, outputfilename):
    """
    explore kaiju protein hits to explore if the protein hits repartition can help to confirm the 
    presence/asbence of a species 
    """
    dict_prot_id={}
    
    for index, row in df_protein.iterrows():
        if "," in row["accession_num"]:
            prot_id = row["accession_num"].split(',')[0]
        else:
            prot_id = row["accession_num"]
        if prot_id in dict_prot_id.keys():
            value = dict_prot_id.get(prot_id)
            value+=1
            dict_prot_id.update({prot_id:value})
        else:
             dict_prot_id.update({prot_id:1})

    msg_prot = "NB of different protein hit on " + filter_term + ": " + str(len(dict_prot_id)) + \
        ", please check proteins.txt file"    
    print(msg_prot)
    
    prot_outputfilename = outputfilename + "_proteins.txt"
    prot_outputfilename2 = outputfilename + "_proteins_full.txt"
    with open(prot_outputfilename, "w") as file:
        file.write(json.dumps(dict_prot_id)) 
   

    df_protein.to_csv(prot_outputfilename2, sep='\t', index=False)


                
            
def filter_taxonomy(df_C_data, ncbi, filter_term, reads_tot, outputfilename, explore_proteins):
    """
    select reads classified according to taxonomy
    """
    dict_read={}
    df_protein = pd.DataFrame()
    list_selected_reads=[]
    invalid_reads = ""
    count_invalid = 0
    count=0
    df_C_data.reset_index(inplace = True) 

    read_name_list = df_C_data["Read_name"].tolist()
    tax_id = df_C_data["ncbi_id"].tolist()
    for i in range(0, len(tax_id)): 
        id = tax_id[i]
        try:
            lineage = ncbi.get_lineage(id)
        except ValueError:
            invalid_reads += "Read: " + str(read_name_list[i]) + " removed because invalid taxid: " + str(id)
            count_invalid += 1
        names = ncbi.get_taxid_translator(lineage)
        list_name=[]
        for (key, value) in names.items():
            list_name.append(value)
        if filter_term in list_name:
            list_selected_reads.append(read_name_list[i])
            count+=1
            
            df_protein = df_protein.append(df_C_data.loc[[i]])
        #dict_read.update( {read_name_list[i]:list_name})
        
        #'HISEQ:105:C4085ANXX:1:2113:4205:42974': ['root', 'Eukaryota', 'Rhodophyta', 'Florideophyceae', 'cellular organisms', 'Corallinophycidae', 'Corallinophycidae sp. GPJ-2013', 'unclassified Corallinophycidae']
    
    if explore_proteins:
        explore_protein_data(df_protein, filter_term, outputfilename)

    outfile = outputfilename + "_Invalid_read_taxid.txt"
    msg_invalid = "NB of reads with invalid taxid " + str(count_invalid) \
        + " " + str(round((count_invalid/reads_tot)*100, 2)) + "%"
    print(msg_invalid)   
    f = open(outfile, "w", encoding="utf-8")
    f.write(invalid_reads) 
    f.close()

    #for (key, value) in dict_read.items():
    #    if filter_term in value:
    #        list_selected_reads.append(key)
    #        count+=1
    msg_filter = "NB of reads with " + filter_term + " classification: " + str(count)\
         + " " + str(round((count/reads_tot)*100, 2)) + "% of all"
    print(msg_filter)  
    msg_filter2 = "NB of reads with " + filter_term + " classification: " + str(count)\
         + " " + str(round((count/len(df_C_data))*100, 2)) + "% of classified"
    print(msg_filter2)  
    
    return list_selected_reads, dict_read

def get_flagged_reads(list_selected_reads, read_file_path1, read_file_path2, paired, filter_term, reads_tot, outputfilename):
    """
    create a read file (fastq) with reads filtered from taxonomic assignation
    """
    index_read_content = []
    if len(list_selected_reads)==0:
        print("No reads in the selection")
        return None
    tmp_list_selected_reads = list_selected_reads
    count=0
    if " " in filter_term:
        filter_term_OUT = filter_term.replace(" ", "_")
    else:
        filter_term_OUT = filter_term
    if paired:
        outfile = outputfilename + "_" + filter_term_OUT + "_R1.fastq"
    else:
        outfile = outputfilename + "_" + filter_term_OUT + ".fastq"
    fi = open(outfile, "w", encoding="utf-8")

    with open(read_file_path1) as f:
        read_data = f.readlines()
        for i in range(0, len(read_data),4): 
            # test optimisation gain en parcourt 1 ligne sur 4
            # Mais perte en faisant readlines ?
            # enumerate ne charge pas en memoire meilleure strategie ?
            line = read_data[i]
        #for (num,line) in enumerate(f):
        #    if line.startswith("@"):
            read_head = str(line.split(" ")[0][1:])
            if read_head in tmp_list_selected_reads: #### list best structure ?
                tmp_list_selected_reads.remove(read_head) # test optimisation
                fi.write(read_data[i])
                fi.write(read_data[i+1])
                fi.write(read_data[i+2])
                fi.write(read_data[i+3])    
                count+=1
                index_read_content.append(i)
                if len(tmp_list_selected_reads) == 0:# test optimisation
                    break      
    msg_readfile = "NB of reads with " + filter_term + " classification in " + outfile \
    + ": " + str(count) + " " + str(round((count/reads_tot)*100, 2)) + "%"
    print(msg_readfile)         
    fi.close()
    
    if paired: ## TODO utiliser index_read_content
        outfile = outputfilename + "_" + filter_term_OUT + "_R2.fastq"
        fi = open(outfile, "w", encoding="utf-8")      
        count=0
        with open(read_file_path2) as f:
            read_data = f.readlines()
            for i in range(0, len(read_data),4): # NB de reads  
            #for i, line in enumerate(read_data): # NB de reads  
                    #if i in index_read_content: # NB de reads to select
                    if count == len(index_read_content):
                        break
                    if i == index_read_content[count]: 
                        # fonctionne QUE si les reads sont ordones pareil dans les deux fichiers (et meme nombre)
                        count+=1
                        fi.write(read_data[i])
                        fi.write(read_data[i+1])
                        fi.write(read_data[i+2])
                        fi.write(read_data[i+3])
                            
                    
        msg_readfile = "NB of reads with " + filter_term + " classification in " + outfile \
        + ": " + str(count) + " " + str(round((count/reads_tot)*100, 2)) + "%"
        print(msg_readfile)   
        fi.close()

if __name__ == "__main__":
    
    ncbi = NCBITaxa()
    #result_file_path = "kaiju.out"
    #read_file_path1 = "/home/jrollin/Desktop/Performance_testing/Sample3/reads_D3.R1.fastq"
    #read_file_path2 = "/home/jrollin/Desktop/Performance_testing/Sample3/reads_D3.R2.fastq"
    #paired = True
    #filter_term = "Viruses"

    parser = argparse.ArgumentParser()
    parser.add_argument('-resultfile', help="output file of the kaiju run", required=True)
    parser.add_argument('-paired',  default=False, action='store_true',\
        help="Are the read in two paired file ?", required=False)
    parser.add_argument('-prot_analysis',  default=False, action='store_true',\
        help="do you want to knonw more about the proteins hits", required=False)
    parser.add_argument('-readfile', help="File containing reads", required=True)
    parser.add_argument('-filepairedread', help="File containing second set of paired reads", required=False)
    parser.add_argument('-taxoclass', help="Taxonomic class hrom which you want to isolate reads", required=True)
    parser.add_argument('-outputfile', help="outputfile ID", required=True)
    args = parser.parse_args()
    result_file_path = args.resultfile
    paired = args.paired
    read_file_path1 = args.readfile
    if paired:
        read_file_path2 = args.filepairedread
    else:
        read_file_path2 = ""
    filter_term = args.taxoclass
    outputfilename = args.outputfile
    prot_analysis = args.prot_analysis

    df_C_data, df_U_data, reads_tot, prot_analysis_possible = get_kaiju_data(result_file_path)
    
    if prot_analysis :
        if prot_analysis_possible:
            explore_proteins=True
        else:
            explore_proteins=False
            print("Proteins exploration of the kaiju output impossible... ")
    else:
        explore_proteins=False

    #print(df_C_data)
    #print(df_U_data)
    # TODO faire quelque chose des U_data => mettre dans un fichiers pour travailler dessus ?
    # cela prendra du temps ...
    list_selected_reads, dict_read = filter_taxonomy(df_C_data, ncbi, filter_term, \
        reads_tot, outputfilename, explore_proteins)
    get_flagged_reads(list_selected_reads, read_file_path1, read_file_path2, paired, \
        filter_term, reads_tot, outputfilename)
