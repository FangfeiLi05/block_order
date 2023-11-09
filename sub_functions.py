import numpy as np
import pandas as pd
import itertools
import csv
import pickle
import psycopg2
import shutil
import os

import Twist_APItools as Twist
import IDT_APItools as IDT



def sliding_window(element, window_size):
    """
    GENERATE A SLIDING WINDOW LIST
    --------------------------------------------------
    INPUTS:
    -- element [str or list]
    -- window_size [int]
    
    OUTPUTS:
    -- [list]: a list of sub-strings or a list of sub-lists
    """
    if len(element) <= window_size:
        return [element]
    else:
        return [element[i:i+window_size] for i in range(len(element) - window_size + 1)]



def convert_blast_pickle_to_txt(filename_pickle, filename_txt):
    """
    CONVERT BLAST RESULT PICKLE FILE TO TXT FILE
    """
    #load pickle file
    with open(filename_pickle, 'rb') as f:
        pickle_file = pickle.load(f)
    f.close()
            
    #write to txt file
    with open(filename_txt, 'w') as f:
        for k in range(len(pickle_file)):
            f.write('\n')
            f.write('score = {}\n'.format(str(pickle_file[k][0])))
            f.write('length = {}\n'.format(str(pickle_file[k][1])))
            f.write('query position = {}-{}\n'.format(str(pickle_file[k][2]), str(pickle_file[k][3])))
            f.write('target position = {}-{}\n'.format(str(pickle_file[k][4]), str(pickle_file[k][5])))
            f.write('percent identity = {}\n'.format(str(pickle_file[k][10])))
            f.write('{}\n'.format(pickle_file[k][7]))
            f.write('{}\n'.format(pickle_file[k][8]))
            f.write('{}\n'.format(pickle_file[k][9]))
    f.close()




def write_to_fasta_or_txt(element, filename):
    """
    WRITE A STRING OR A LIST OF STRING TO A FILE (txt or fasta)
    --------------------------------------------------
    INPUTS:
    -- element [str or list]
    -- filename [str]
        
    OUTPUTS:
    --
    """
    if type(element) == str:
        element = [element]
            
    with open(filename, 'w') as f:
        for k in range(len(element)):
            f.write('>{}\n'.format(k))
            f.write('{}\n'.format(element[k]))
    f.close()



def write_dict_to_csv(data_dict, filename, option):
    """
    WRITE DICT TO A FILE (csv or xlsx)
    --------------------------------------------------
    INPUTS:
    -- data_dict [dict]
    -- filename [str]
    -- option [str]: 'csv', 'xlsx'
    OUTPUTS:
    --
    """
    tmp = list(itertools.zip_longest(*list(data_dict.values())))
    if option == 'csv':
        with open(filename, 'w') as f:
            w = csv.writer(f)
            w.writerow(data_dict.keys())
            w.writerows(tmp)
        f.close()
    
    elif option == 'xlsx':
        df = pd.DataFrame(tmp)
        with pd.ExcelWriter(filename) as f:
            df.to_excel(f, header=list(data_dict.keys()), index = False)

    #df = pd.DataFrame(data_dict)
    #with pd.ExcelWriter('{}_order_IDT.xlsx'.format(filename_pre)) as f:
    #    df.to_excel(f, index = False)


#worked, but not yet carefully read line by line
def twist_buildable(sequence_list, label_list, output_analysis=False, filename_pre=''):
    """
    CHECK WHETHER TWIST CAN BUILD
    --------------------------------------------------
    INPUTS:
    -- sequence_list [list]: a list of sequence to be checked
    -- label_list [list]: a list of label corresponding for each of the sequence
    
    OUTPUTS:
    -- twist_result [list]: a list of lists [[1], [2], [3]]:
    1. TRUE/FALSE: whether Twist can build the construct or not
    2. Twist difficulty rating
    3. Twist build  issues identified
    """
    len_sequence_list = len(sequence_list)
    constructs_result = [0] * len_sequence_list
    twist_score_result = [0] * len_sequence_list
    twist_issue_result = [0] * len_sequence_list
        
    for k in range(len_sequence_list):
        constructs_result[k] = Twist.submit_gene_block(sequence_list[k], label_list[k])
        #batch all sumbissions to avoid waiting for each score
                
        twist_score_result[k] = Twist.get_construct_score(constructs_result[k])
        #get all scoring info
    
    ######
    if output_analysis:
        with open('{}_twist_analysis.pickle'.format(filename_pre), 'wb') as f:
            pickle.dump(twist_score_result, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.close()
    ######
    
    twist_buildable_result = [i['score'] == 'BUILDABLE' for i in twist_score_result]
    #determine which blocks are buildable
    twist_difficulty_result = [i['score_data']['difficulty'] for i in twist_score_result]
    #determine the score for each block
        
    for i in range(len_sequence_list):
        if len(twist_score_result[i]['score_data']['issues']) == 0: #without issue (check whether a list is empty, i.e., [])
            twist_issue_result[i] = 'None'
        else: #with issue
            tmp_issue_1 = twist_score_result[i]['score_data']['issues']
            tmp_issue_2 = []
            for j in range(len(tmp_issue_1)):
                if tmp_issue_1[j]['storedAt'] == None:
                    tmp_issue_2.append('{} {}'.format(tmp_issue_1[j]['code'],
                                                          tmp_issue_1[j]['title']))
                else:
                    tmp_issue_2.append('{} {} {}-{}'.format(tmp_issue_1[j]['code'],
                                                                tmp_issue_1[j]['title'],
                                                                tmp_issue_1[j]['location']['start'],
                                                                tmp_issue_1[j]['location']['end']))
                twist_issue_result[i] = '; '.join(tmp_issue_2)
        
    twist_result = [twist_buildable_result, twist_difficulty_result, twist_issue_result]
    
    return twist_result





def idt_buildable(sequence_list, label_list, output_analysis=False, filename_pre=''):
    """
    CHECK WHETHER IDT CAN BUILD
    --------------------------------------------------
    INPUTS:
    -- sequence_list [list]: a list of sequence to be checked
    -- label_list [list]: a list of label corresponding for each of the sequence
    
    OUTPUTS:
    -- idt_result [list]: a list of lists [[1], [2], [3]]:
    1. TRUE/FALSE: whether IDT can build the construct or not
    2. IDT difficulty rating (higher means worse)
    3. IDT build issues identified
    """
    
    token = IDT.get_access_token("***", "***", "***", "***") #Get access token (lasts ~1h)
    len_sequence_list = len(sequence_list)
    idt_analysis_result = [0] * len_sequence_list
    idt_scores_result = [0] * len_sequence_list
    idt_buildable_result = [0] * len_sequence_list
    idt_issues_result = [0] * len_sequence_list
        
    for i in range(len_sequence_list):
        idt_analysis_result[i] = IDT.Gblock_complexity(sequence_list[i], token)
        
    if output_analysis:
        with open('{}_idt_analysis.pickle'.format(filename_pre), 'wb') as f:
            pickle.dump(idt_analysis_result, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.close()
               
    for i in range(len_sequence_list):
        tmp_analysis = idt_analysis_result[i]
        if type(tmp_analysis) == dict:
            idt_scores_result[i] = tmp_analysis['Score']
            idt_buildable_result[i] = tmp_analysis['Score'] < 10 #!!!why 10
            idt_issues_result[i] = tmp_analysis['DisplayText']
        elif type(tmp_analysis) == list:
            tmp = sum([j['Score'] for j in tmp_analysis])
            idt_scores_result[i] = tmp
            idt_buildable_result[i] = tmp < 10
            idt_issues_result[i] = '; '.join([j['DisplayText'] for j in tmp_analysis])
        
    idt_result = [idt_buildable_result, idt_scores_result, idt_issues_result]
                
    return idt_result
        




def output_to_worksheet(block_file, filename_pre):
    """
    INPUTS:
    block_file [str]: a .csv file
    
    OUTPUTS:
    1. three xlsx files for ordering from Twist, IDT and Genscript
    2.
    """
    csv_input = pd.read_csv(block_file)
    block_dict = {ele: list(csv_input[ele]) for ele in list(csv_input.keys())}
    block_num = len(csv_input)
    
    idx_twist = [True if block_dict['Twist_buildable'][k] else False for k in range(block_num)]
    
    idx_idt = [True if ((not block_dict['Twist_buildable'][k])
                    and (block_dict['IDT_buildable'][k]))
               else False for k in range(block_num)]
               
    idx_genscript = [True if ((not block_dict['Twist_buildable'][k])
                          and (not block_dict['IDT_buildable'][k]))
                     else False for k in range(block_num)]
                     
    block_label_array = np.array([block_dict['DNA_block_name'][k] for k in range(block_num)])
    block_sequence_array = np.array([block_dict['DNA_block_sequence'][k] for k in range(block_num)])
    
    supplier_array = np.zeros(block_num).astype('<U100')
    supplier_array[idx_twist] = 'Twist'
    supplier_array[idx_idt] = 'IDT'
    supplier_array[idx_genscript] = 'GenScript'
    supplier_list = list(supplier_array)
    
    ##### save the blockOrder csv file
    blockOrder_database_dict = {'alias': block_dict['DNA_block_name'],
                                'seq': block_dict['block_without_adapter'],
                                'adapterLeft': block_dict['adapter_left'],
                                'adapterRight': block_dict['adapter_right'],
                                'hash': block_dict['DNA_block_sequence_hash'],
                                'supplier': supplier_list,
                                'person': [],
                                'purchasedAt': [],
                                'receivedAt': [],
                                'concentration': [],
                                'volumn': [],
                                'storedAt': []}
                                 
    filename_blockOrder_database_csv = '{}_blockOrder_database.csv'.format(filename_pre)
    write_dict_to_csv(blockOrder_database_dict, filename_blockOrder_database_csv, 'csv')
    
    #####
    filename_blockOrder_database_csv_tmp = '/tmp/{}_blockOrder_database.csv'.format(filename_pre.split('/')[-1])
    shutil.move(filename_blockOrder_database_csv, filename_blockOrder_database_csv_tmp)
    
    change_right = 'chmod 777 {}'.format(filename_blockOrder_database_csv_tmp)
    os.system(change_right)

    conn = psycopg2.connect(database='my_database',
                            user='postgres',
                            password='870525',
                            host='localhost',
                            port='5432')
                            
    #conn = psycopg2.connect(database='our_database',
    #                        user='postgres',
    #                        password='pro3165',
    #                        host='localhost',
    #                        port='5432')

    cur = conn.cursor()
    cur.execute('CALL queryfromblockOrder (%s)', (filename_blockOrder_database_csv_tmp,))
    cur.close()
    conn.close()
    
    shutil.move(filename_blockOrder_database_csv_tmp, filename_blockOrder_database_csv)

    
    csv_input = pd.read_csv(filename_blockOrder_database_csv)
    blockOrder_dict = {ele: list(csv_input[ele]) for ele in list(csv_input.keys())}
    #location_list = blockOrder_dict['storedAt']
    location_list = blockOrder_dict['storedat']
    
    idx_in_database = [type(location_list[k])==str for k in range(block_num)]
    idx_notin_database = [not type(location_list[k])==str for k in range(block_num)]
    
    if sum(idx_twist):
        twist_dict = {'Sequence name':block_label_array[idx_twist],
                      'Sequence':block_sequence_array[idx_twist]}
        df = pd.DataFrame(twist_dict)
        with pd.ExcelWriter('{}_order_Twist.xlsx'.format(filename_pre)) as f:
            df.to_excel(f, index = False)
        
    if sum(idx_idt):
        idx_dict = {'Name':block_label_array[idx_idt],
                    'Sequence':block_sequence_array[idx_idt],
                    '5p Phosphorylation':[]}
        filename_xlsx_idx = '{}_order_IDT.xlsx'.format(filename_pre)
        write_dict_to_csv(idx_dict, filename_xlsx_idx, 'xlsx')
        
    if sum(idx_genscript):
        genscript_dict = {'Gene name':block_label_array[idx_genscript],
                          'Gene sequence':block_sequence_array[idx_genscript],
                          '5p Flanking Sequence':[],
                          '3p Flanking Sequence':[]}
        filename_xlsx_genscript = '{}_order_Genscript.xlsx'.format(filename_pre)
        write_dict_to_csv(genscript_dict, filename_xlsx_genscript, 'xlsx')
            
    if sum(idx_in_database) > 0:
        idx_twist_notin_database = [idx_twist[k] and idx_notin_database[k] for k in range(block_num)]
        if sum(idx_twist_notin_database):
            twist_dict = {'Sequence name':block_label_array[idx_twist_notin_database],
                          'Sequence':block_sequence_array[idx_twist_notin_database]}
            df = pd.DataFrame(twist_dict)
            with pd.ExcelWriter('{}_order_Twist_not_in_database.xlsx'.format(filename_pre)) as f:
                df.to_excel(f, index = False)
        
        idx_idt_notin_database = [idx_idt[k] and idx_notin_database[k] for k in range(block_num)]
        if sum(idx_idt_notin_database):
            idx_dict = {'Name':block_label_array[idx_idt_notin_database],
                        'Sequence':block_sequence_array[idx_idt_notin_database],
                        '5p Phosphorylation':[]}
            filename_xlsx_idx = '{}_order_IDT_not_in_database.xlsx'.format(filename_pre)
            write_dict_to_csv(idx_dict, filename_xlsx_idx, 'xlsx')
            
        idx_genscript_notin_database = [idx_genscript[k] and idx_notin_database[k] for k in range(block_num)]
        if sum(idx_genscript_notin_database):
            genscript_dict = {'Gene name':block_label_array[idx_genscript_notin_database],
                              'Gene sequence':block_sequence_array[idx_genscript_notin_database],
                              '5p Flanking Sequence':[],
                              '3p Flanking Sequence':[]}
            filename_xlsx_genscript = '{}_order_Genscript_not_in_database.xlsx'.format(filename_pre)
            write_dict_to_csv(genscript_dict, filename_xlsx_genscript, 'xlsx')
       
    

    cloning_worksheet_dict = {'Experiment_ID': block_dict['experiment_id'],
                              'DNA_block_name': block_dict['DNA_block_name'],
                              'Supplier': supplier_list,
                              'Plasmid_backbone': block_dict['donor_plasmid'],
                              'Marker':  block_dict['donor_plasmid_marker'],
                              'Marker_counter':  block_dict['donor_plasmid_marker_counter'],
                              'Built (True/False)':[],
                              'Success_rate':[],
                              'Turn_around_time (days)':[],
                              'Gap_repair (True/False)':[],
                              'Cut/ligate (True/False)':[],
                              'Gibson_assembly (True/False)':[]}
    
    filename_csv = '{}_cloning_worksheet.csv'.format(filename_pre)
    write_dict_to_csv(cloning_worksheet_dict, filename_csv, 'csv')
    
    #filename_xlsx = '{}_cloning_worksheet.xlsx'.format(filename_pre)
    #write_dict_to_csv(cloning_worksheet_dict, filename_xlsx, 'xlsx')
    
    
    #####
    
    
    #####
    
    
def my_print(str, option):
    if option:
        print(str)
