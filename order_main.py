import numpy as np
import pandas as pd
import csv
import os
import json
from Bio import SeqIO
from multiprocess import Pool
#from tqdm import tqdm
import time
import argparse

import find_cut
import calculate_score
import return_block
import process_file
from sub_functions import my_print

import matplotlib.pyplot as plt




class OrderMain(object):
    """
    """
    def __init__(self, parameter_csv):
        """
        """
        csv_input = pd.read_csv(parameter_csv)
        self.param_dict = {param: value for param, value
                                        in zip(csv_input['parameter'], csv_input['value']) }
        
        self.parallel = json.loads(self.param_dict['parallel'].lower())
        #self.plot = json.loads(self.param_dict['plot'].lower())
        #whether make and save a plot

        self.l_min = int(self.param_dict['min_length'])
        self.l_max = int(self.param_dict['max_length'])
        self.o_min = int(self.param_dict['min_overlapping'])
        self.o_max = int(self.param_dict['max_overlapping'])
        self.n_max = int(self.param_dict['max_number'])
            
        self.dir_output = self.param_dict['output_directory']
        
        score_option = self.param_dict['score_option']
        if 'BLAST' in score_option: #BLAST
            use_BLAST = True
        else:
            use_BLAST = False
            
        if use_BLAST:
            dir_blast = '{}blast_result/'.format(self.dir_output)
            #folder for BLAST result
            dir_blast_database = '{}database/'.format(dir_blast)
            #folder for database in BLAST
            dir_blast_tmp = '{}tmp/'.format(dir_blast)
            #folder for temporary in BLAST
            if not os.path.exists(dir_blast):
                os.mkdir(dir_blast)
            if not os.path.exists(dir_blast_database):
                os.mkdir(dir_blast_database)
            if not os.path.exists(dir_blast_tmp):
                os.mkdir(dir_blast_tmp)
                
                
                
        #self.dir_score = '{}position_score/'.format(self.dir_output)
        #if not os.path.exists(self.dir_score):
        #    os.mkdir(self.dir_score)
        
        #if self.plot:
        #    self.dir_plot = '{}gene_plots/'.format(self.dir_output)
        #    if not os.path.exists(self.dir_plot):
        #       os.mkdir(self.dir_plot)
        
        
    
    def process_target_single(self, sequence, label):
        """
        FUNCTION FOR PROCESSING EACH DNA SEQUENCE IN A FASTA FILE
        ----------------------------------------
        INPUTS:
        -- option [str]: 'all', 'score_only'
        """
        my_print('   -- Calculating position score...', not self.parallel)
        score_dict = calculate_score.score(sequence, self.param_dict, label)
        score_result = score_dict['summed_score']
        
        ####################
        #print('aaa')
        #print( score_dict['primer3_hairpin_score']['norm'])
        
        fig = plt.figure()
        plt.plot(np.arange(len(sequence)), score_dict['gc_score']['norm'], '.', color='red')
        plt.plot(np.arange(len(sequence)), score_dict['blast_score']['norm'], '.', color='blue')
        plt.plot(np.arange(len(sequence)), score_dict['primer3_hairpin_score']['norm'], '.' , color='orange')
        plt.plot(np.arange(len(sequence)), score_dict['ndu_repeat_score']['norm'], '.', color='green')
        
        
        
        my_print('\n   -- Finding optimal cut...', not self.parallel)
        cut_dict = find_cut.cut(score_result, n_max = self.n_max,
                                              l_min = self.l_min,
                                              l_max = self.l_max,
                                              o_min = self.o_min,
                                              o_max = self.o_max,
                                              print_option = not self.parallel)
                                                      
        cut_result = [cut_dict['start_pos'], cut_dict['end_pos']]
        
        my_print('\n   -- Returning output block...', not self.parallel)
        block_result = return_block.block(sequence, cut_result, self.param_dict, label)
        
        return block_result
        

        

def submain(parameter_csv):
    """
    """
    csv_input = pd.read_csv(parameter_csv)
    param_dict = {param: value for param, value
                               in zip(csv_input['parameter'], csv_input['value']) }
    
    experiment_id = param_dict['experiment_id']
    input_file = param_dict['input_file']
    
    dir_output = param_dict['output_directory']
    dir_blast = '{}blast_result/'.format(dir_output)
    dir_blast_tmp = '{}tmp/'.format(dir_blast)
    dir_worksheet = '{}worksheet/'.format(dir_output)
    if not os.path.exists(dir_worksheet):
        os.mkdir(dir_worksheet)
            
    print('Starting for: {} ...'.format(experiment_id))
    t_start = time.time() ###
    f = open('{}{}_DNA_block_info.csv'.format(dir_worksheet, experiment_id), 'w')
    w = csv.writer(f)
    w.writerow(['experiment_id',
                'DNA_block_name',
                'DNA_label',
                'block_number',
                'DNA_block_sequence',
                'DNA_block_sequence_hash',
                'block_length',
                'block_start',
                'block_end',
                'donor_plasmid',
                'donor_plasmid_marker',
                'donor_plasmid_marker_counter',
                'block_without_adapter',
                'adapter_left',
                'adapter_right'])
                
    my_object = OrderMain(parameter_csv)
    for record in SeqIO.parse(input_file, 'fasta'):
        label = record.id
        sequence = str(record.seq)
        print('\n----------------------------------------' )
        print('-- For sequence with label "{}"...'.format(label))
        output_single = my_object.process_target_single(sequence, label)
        w.writerows(output_single)
    f.close()
    t_end = time.time() ###
    print('Total time for 1st part: {}'.format(t_end - t_start)) ###
    
    
    print('\n----------------------------------------' )
    print('Writing to worksheet...')
    t_start = time.time() ###
    block_filename = '{}{}_DNA_block_info.csv'.format(dir_worksheet, experiment_id)
    process_file.file(block_filename, param_dict)
    
    os.system('rm -r -f {}'.format(dir_blast_tmp))
    print('Done!' )
    t_end = time.time() ###
    print('Total time for 2nd part: {}'.format(t_end - t_start)) ###
    
    


def submain_parallel(parameter_csv):
    """
    """
    csv_input = pd.read_csv(parameter_csv)
    param_dict = {param: value for param, value
                               in zip(csv_input['parameter'], csv_input['value']) }
    
    experiment_id = param_dict['experiment_id']
    input_file = param_dict['input_file']
        
    dir_output = param_dict['output_directory']
    dir_blast = '{}blast_result/'.format(dir_output)
    dir_blast_tmp = '{}tmp/'.format(dir_blast)
    dir_worksheet = '{}worksheet/'.format(dir_output)
    if not os.path.exists(dir_worksheet):
        os.mkdir(dir_worksheet)
            
    print('Starting for: {} ...'.format(experiment_id))
    t_start = time.time() ###
    label_list = []
    sequence_list = []
    for record in SeqIO.parse(input_file, 'fasta'):
        label_list.append(record.id)
        sequence_list.append(str(record.seq))
        
    ##########
    my_object = OrderMain(parameter_csv)
    my_obj_pool = Pool()  #Pool(processes=n), default: number of logical CPU cores, print(my_obj_pool._processes)
    args = [(sequence_list[k], label_list[k]) for k in range(len(label_list))]
    output_all = my_obj_pool.starmap(my_object.process_target_single, args)
    my_obj_pool.close()
    
    #my_obj_pool = Pool(processes=8)
    #output_async = [my_obj_pool.apply_async(my_object.process_target_single,
    #                                        args=(sequence_list[k],
    #                                              label_list[k]))
    #                for k in range(len(label_list))]
    #output_all = [ar.get() for ar in output_async]
    #my_obj_pool.close()
    
    f = open('{}{}_DNA_block_info.csv'.format(dir_worksheet, experiment_id), 'w')
    w = csv.writer(f)
    w.writerow(['experiment_id',
                'DNA_block_name',
                'DNA_label',
                'block_number',
                'DNA_block_sequence',
                'DNA_block_sequence_hash',
                'block_length',
                'block_start',
                'block_end',
                'donor_plasmid',
                'donor_plasmid_marker',
                'donor_plasmid_marker_counter',
                'block_without_adapter',
                'adapter_left',
                'adapter_right'])
    for output_single in output_all:
        w.writerows(output_single)
    f.close()
    t_end = time.time() ###
    print('Total time for 1st part: {}'.format(t_end - t_start)) ###
    
    print('\n----------------------------------------' )
    print('Writing to worksheet...')
    t_start = time.time() ###
    block_filename = '{}{}_DNA_block_info.csv'.format(dir_worksheet, experiment_id)
    process_file.file(block_filename, param_dict)
    
    #os.system('rm -r -f {}'.format(dir_blast_tmp))
    print('Done!' )
    t_end = time.time() ###
    print('Total time for 2nd part: {}'.format(t_end - t_start)) ###


def process(parameter_csv):
    """
    """
    csv_input = pd.read_csv(parameter_csv)
    param_dict = {param: value for param, value
                               in zip(csv_input['parameter'], csv_input['value']) }
    
    parallel = json.loads(param_dict['parallel'].lower())

    if parallel:
        submain_parallel(parameter_csv)
    elif not parallel:
        submain(parameter_csv)
        
    
    


def main():
    """
    """
    parser = argparse.ArgumentParser(description='Process DNA sequence',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-p', '--parameter', type=str, required=True,
                        help='a .csv file with parameters')
                    
    args = parser.parse_args()
    
    parameter_csv = args.parameter
    csv_input = pd.read_csv(parameter_csv)
    param_dict = {param: value for param, value
                               in zip(csv_input['parameter'], csv_input['value']) }
    
    parallel = json.loads(param_dict['parallel'].lower())

    if parallel:
        submain_parallel(parameter_csv)
    elif not parallel:
        submain(parameter_csv)
    
    
if __name__ == "__main__":
    main()

    
