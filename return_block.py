import numpy as np
import os
import yaml
import hashlib


class ReturnBlock(object):
    def __init__(self, sequence, cut_pos, param_dict):
        
        self.experiment_id = param_dict['experiment_id']
        #the assembly batch, each batch must have a unique ID
        
        #donor plasmid
        self.first_donor_plasmid = param_dict['first_donor_plasmid']
        #1st odd donor plasmid (used only once)
        self.even_donor_plasmid = param_dict['even_donor_plasmid']
        #even donor plasmid (used every 2 stitches)
        self.odd_donor_plasmid = param_dict['odd_donor_plasmid']
        #odd donor plasmid (used every 2 stitches)
        
        #selectable markers on each plasmid
        file_yaml_donor = param_dict['yaml_donor']
        yaml_donor_input = open(file_yaml_donor, 'r')
        param_donor_dict = yaml.load(yaml_donor_input, Loader=yaml.FullLoader)
        self.first_donor_plasmid_marker = param_donor_dict[self.first_donor_plasmid]['marker']
        self.first_donor_plasmid_marker_counter = param_donor_dict[self.first_donor_plasmid]['marker_counter']
        self.first_donor_header = param_donor_dict[self.first_donor_plasmid]['header_seq']
        self.first_donor_tail = param_donor_dict[self.first_donor_plasmid]['tail_seq']
        
        self.even_donor_plasmid_marker = param_donor_dict[self.even_donor_plasmid]['marker']
        self.even_donor_plasmid_marker_counter = param_donor_dict[self.even_donor_plasmid]['marker_counter']
        self.even_donor_header = param_donor_dict[self.even_donor_plasmid]['header_seq']
        self.even_donor_tail = param_donor_dict[self.even_donor_plasmid]['tail_seq']
        
        self.odd_donor_plasmid_marker = param_donor_dict[self.odd_donor_plasmid]['marker']
        self.odd_donor_plasmid_marker_counter = param_donor_dict[self.odd_donor_plasmid]['marker_counter']
        self.odd_donor_header = param_donor_dict[self.odd_donor_plasmid]['header_seq']
        self.odd_donor_tail = param_donor_dict[self.odd_donor_plasmid]['tail_seq']
        

        self.gene_tail_extension = param_dict['gene_tail_extension']
        #extra sequence to place at the end in case a buffer is needed by the company
                
                
    def process_result(self, sequence, cut_pos, label):
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
        
        sequence_pad = sequence + self.gene_tail_extension
        pos_start, pos_end = cut_pos[0], cut_pos[1]
        len_pos = len(pos_start)
        
        block_label_list = []
        block_sequence_list = []
        block_donor_plasmid_list = []
        block_donor_plasmid_marker_list = []
        block_donor_plasmid_marker_counter_list = []
        adapter_left_list = []
        adapter_right_list = []
        block_without_adapter_list = []
        
        for k in range(len_pos):
            block_label_list.append(label[0:25] + '_' + str(k+1)) #needs to be short enough for twist API to accept it,

            tmp = sequence_pad[pos_start[k]-1:pos_end[k]]
            block_without_adapter_list.append(tmp)
                            
            if k == 0: #first odd_donor
                adapter_left_list.append(self.first_donor_header)
                adapter_right_list.append(self.first_donor_tail)
                block_sequence_list.append(self.first_donor_header + tmp + self.first_donor_tail)
                block_donor_plasmid_list.append(self.first_donor_plasmid)
                block_donor_plasmid_marker_list.append(self.first_donor_plasmid_marker)
                block_donor_plasmid_marker_counter_list.append(self.first_donor_plasmid_marker_counter)
                
            elif np.mod(k, 2) == 1: #even_donor
                adapter_left_list.append(self.even_donor_header)
                adapter_right_list.append(self.even_donor_tail)
                block_sequence_list.append(self.even_donor_header + tmp + self.even_donor_tail)
                block_donor_plasmid_list.append(self.even_donor_plasmid)
                block_donor_plasmid_marker_list.append(self.even_donor_plasmid_marker)
                block_donor_plasmid_marker_counter_list.append(self.even_donor_plasmid_marker_counter)
                
            elif np.mod(k, 2) == 0: #odd_donor
                adapter_left_list.append(self.odd_donor_header)
                adapter_right_list.append(self.odd_donor_tail)
                block_sequence_list.append(self.odd_donor_header + tmp + self.odd_donor_tail)
                block_donor_plasmid_list.append(self.odd_donor_plasmid)
                block_donor_plasmid_marker_list.append(self.odd_donor_plasmid_marker)
                block_donor_plasmid_marker_counter_list.append(self.odd_donor_plasmid_marker_counter)
        
        
        block_sequence_list = [s.upper() for s in block_sequence_list]
        block_without_adapter_list = [s.upper() for s in block_without_adapter_list]
        adapter_left_list = [s.upper() for s in adapter_left_list]
        adapter_right_list = [s.upper() for s in adapter_right_list]
                                         
        block_result = []
        for k in range(len_pos):
            h = hashlib.sha256()
            h.update(block_sequence_list[k].encode())
            hash_str = h.hexdigest()
            block_result.append([self.experiment_id,
                                 block_label_list[k],
                                 label,
                                 k+1,
                                 block_sequence_list[k],
                                 hash_str,
                                 len(block_sequence_list[k]),
                                 pos_start[k],
                                 pos_end[k],
                                 block_donor_plasmid_list[k],
                                 block_donor_plasmid_marker_list[k],
                                 block_donor_plasmid_marker_counter_list[k],
                                 block_without_adapter_list[k],
                                 adapter_left_list[k],
                                 adapter_right_list[k]]) #might need to add more
        return block_result
        
        
        

def block(sequence, cut_pos, param_dict, label):
    """
    """
    my_object = ReturnBlock(sequence, cut_pos, param_dict)
    result = my_object.process_result(sequence, cut_pos, label)
    
    return result
