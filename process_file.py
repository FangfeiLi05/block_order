import numpy as np
import pandas as pd
import os
#import json
#from multiprocess import Pool

from sub_functions import twist_buildable, idt_buildable, output_to_worksheet, write_dict_to_csv



class ProcessFile(object):
    def __init__(self, param_dict):
        
        self.experiment_id = param_dict['experiment_id']
        #the assembly batch, each batch must have a unique ID
        
        #self.output_twist_analysis = json.loads(param_dict['Twist_analysis'].lower())
        #self.output_idt_analysis = json.loads(param_dict['IDT_analysis'].lower())
        self.output_twist_analysis = True
        self.output_idt_analysis = True
        
        self.dir_output = param_dict['output_directory']
        
        self.dir_worksheet = '{}worksheet/'.format(self.dir_output)
        if not os.path.exists(self.dir_worksheet):
            os.mkdir(self.dir_worksheet)
        
        if self.output_twist_analysis or self.output_idt_analysis:
            self.dir_vendor = '{}vendor_analysis/'.format(self.dir_output)
            if not os.path.exists(self.dir_vendor):
                os.mkdir(self.dir_vendor)
        
        #self.parallel = json.loads(param_dict['parallel'].lower())
    
    
    
    def process_file(self, block_filename):
        """
        INPUTS:
        -- sequence [str]:
        -- cut_pos [list] : A list of two lists [[1], [2]]:
           1. start position for each block (index starts as 1,2,...)
           2. end position for each block
       
        OUTPUTS:
        A list of lists. Each list reports:
        1) The unique DNA block name
        2) The construct name
        3) The build position for a DNA block. 1 = the first block in an assembly
        4) The sequence of the DNA block
        5) Can Twist build it?
        6) Can IDT build it?
        """
        csv_input = pd.read_csv(block_filename)
        block_dict = {ele: list(csv_input[ele]) for ele in list(csv_input.keys())}
        block_label_list  = block_dict['DNA_block_name']
        block_sequence_list= block_dict['DNA_block_sequence']
        len_list = len(block_label_list)
        
        if self.output_twist_analysis:
            print('-- Determining if Twist can build each block...')
            filename_pre = '{}{}'.format(self.dir_vendor, self.experiment_id)
        
            twist_result = twist_buildable(block_sequence_list,
                                           block_label_list,
                                           self.output_twist_analysis,
                                           filename_pre)
           
            block_dict['Twist_buildable'] = [twist_result[0][k] for k in range(len_list)]
            block_dict['Twist_difficulty'] = [twist_result[1][k] for k in range(len_list)]
            block_dict['Twist_issues'] = [twist_result[2][k] for k in range(len_list)]
        
        else:
            block_dict['Twist_buildable'] = [False for k in range(len_list)]
        
        if self.output_idt_analysis:
            print('-- Determining if IDT can build each block...')
            idt_result = idt_buildable(block_sequence_list,
                                       block_label_list,
                                       self.output_idt_analysis,
                                       filename_pre)
                                       
            block_dict['IDT_buildable'] = [idt_result[0][k] for k in range(len_list)]
            block_dict['IDT_score'] = [idt_result[1][k] for k in range(len_list)]
            block_dict['IDT_issues'] = [idt_result[2][k] for k in range(len_list)]
        
        else:
            block_dict['IDT_buildable'] = [False for k in range(len_list)]
        
            #my_obj_pool = Pool()
            #args = [(block_sequence_list[k],
            #         block_label_list[k],
            #         self.output_idt_analysis,
            #         filename_pre) for k in range(len(block_label_list))]
            #output_all = my_obj_pool.starmap(idt_buildable, args)
            #my_obj_pool.close()
        
        write_dict_to_csv(block_dict, block_filename, 'csv')
        
        filename_pre = '{}{}'.format(self.dir_worksheet, self.experiment_id)
        output_to_worksheet(block_filename, filename_pre) #cloning worksheet, and database
        
        

def file(block_filename, param_dict):
    """
    """
    my_object = ProcessFile(param_dict)
    my_object.process_file(block_filename)

