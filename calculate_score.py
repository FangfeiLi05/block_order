import numpy as np
import yaml
import json
import os
import pickle
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Blast import NCBIXML
import primer3 #pip install primer3-py (https://pypi.org/project/primer3-py/#files)
import DNANDU as DU
from sub_functions import convert_blast_pickle_to_txt, sliding_window, write_to_fasta_or_txt, my_print



class CalculateScore(object):
    def __init__(self, sequence, param_dict, label):
        
        file_yaml_score = param_dict['yaml_score']
        yaml_score_input = open(file_yaml_score, 'r')
        #----------
        param_score_dict = yaml.load(yaml_score_input, Loader=yaml.FullLoader)
        self.GC_window_size = param_score_dict['GC']['GC_window_size']
        self.GC_threshold_upper = param_score_dict['GC']['GC_threshold_upper']
        self.GC_threshold_lower = param_score_dict['GC']['GC_threshold_lower']
        #----------
        self.WORD_SIZE = param_score_dict['BLAST']['WORD_SIZE']
        self.PERC_IDENTITY = param_score_dict['BLAST']['PERC_IDENTITY']
        self.EVALUE = param_score_dict['BLAST']['EVALUE']
        self.REWARD = param_score_dict['BLAST']['REWARD']
        self.PENALTY = param_score_dict['BLAST']['PENALTY']
        self.GAPOPEN = param_score_dict['BLAST']['GAPOPEN']
        self.GAPEXTEND = param_score_dict['BLAST']['GAPEXTEND']
        #----------
        self.hairpin_window_size = param_score_dict['Primer3']['hairpin_window_size']
        self.hairpin_threshold = param_score_dict['Primer3']['hairpin_threshold']
        #----------
        self.NDU_repeat_window_size = param_score_dict['NDU']['NDU_repeat_window_size']
        self.NDU_threshold = param_score_dict['NDU']['NDU_threshold']
        #---------- this option is currently not added yet
        self.MFE_window_size = param_score_dict['MFE']['MFE_window_size']
        self.MFE_skip = param_score_dict['MFE']['MFE_skip']
        self.MFE_offset = param_score_dict['MFE']['MFE_offset']
        self.MFE_threshold = param_score_dict['MFE']['MFE_threshold']

        
        self.recipient_plasmid_ends_file = param_dict['recipient_plasmid_ends_file']
        #file with plasmid sequences surronding assembly blocks (used for DNADU repeat)
        self.plasmids_file = param_dict['plasmids_file']
        #file with all unique sequences in donor and recipient plasmids (used for BLAST)
        
        self.dir_output = param_dict['output_directory']
        #self.dir_score = '{}position_score/'.format(self.dir_output)
        
        #whether use this method to calculate overlap regions (True/False) ....
        score_option = param_dict['score_option']
        if 'GC' in score_option: #GC content
            self.use_GC = True
        else:
            self.use_GC = False
            
        if 'BLAST' in score_option: #BLAST
            self.use_BLAST = True
            self.dir_blast = '{}blast_result/'.format(self.dir_output)
            #folder for BLAST result
            self.dir_blast_database = '{}database/'.format(self.dir_blast)
            #folder for database in BLAST
            self.dir_blast_tmp = '{}tmp/'.format(self.dir_blast)
            #folder for temporary in BLAST
        else:
            self.use_BLAST = False
            
        if 'Primer3' in score_option: #hairpin score in Primer3
            self.use_Primer3 = True
        else:
            self.use_Primer3 = False
            
        if 'NDU' in score_option: #DNADU repeat score (computationally expensive)
            self.use_NDU = True
        else:
            self.use_NDU = False
            
        if 'MFE' in score_option: #minimum free energy (computationally expensive)
            self.use_MFE = True
        else:
            self.use_MFE = False
        
        
        self.label = label
        self.parallel = json.loads(param_dict['parallel'].lower())
        
        
        
    def calcu_gc_content(self, sequence):
        """
        CALCULATE GC-CONTENT CONSIDERING SLIDING WINDOWS
        --------------------------------------------------
        INPUTS:
        -- sequence [str]: a DNA sequence
        
        OUTPUTS:
        -- gc_contents [list]: a list of GC content
        """
        #window_size = self.GC_window_size
        #len_left = int(np.floor(window_size/2))
        #len_right = window_size - len_left - 1
        #sequence_pad = sequence[-len_left:] + sequence + sequence[:len_right]
        #sequence_slide_window = sliding_window(sequence_pad, window_size)
        #gc_score = [gc_fraction(seq)/100 for seq in sequence_slide_window]
        
        window_size = self.GC_window_size
        sequence_slide_window = sliding_window(sequence, window_size)
        tmp = int((window_size - 1)/2)
        gc_content = [0]*tmp + [gc_fraction(seq) for seq in sequence_slide_window] + [0]*tmp
        
        return gc_content
        
        
        
    def make_blast_database(self, sequence_database):
        """
        MAKE A DATABASE
        --------------------------------------------------
        INPUTS:
        -- sequence_database [list]: a list of DNA sequences (target DNA sequence, plasmids (donor and recipient))
        """
        #Make a fatsa file for building BLAST database
        filename = '{}tmp_{}.fa'.format(self.dir_blast_tmp, self.label)
        write_to_fasta_or_txt(sequence_database, filename)
        
        #Make a BLAST database
        os.system('makeblastdb -in {0}tmp_{2}.fa -title {2} -dbtype nucl -input_type fasta -hash_index -out {1}{2} > {1}log_file_{2}.txt'.format(self.dir_blast_tmp, self.dir_blast_database, self.label))
        
        #remove the fatsa file that used for building BLAST database
        os.system('rm -f {}tmp_{}.fa'.format(self.dir_blast_tmp, self.label))
        
    
        
    def make_blast_single(self, query):
        """
        NBLASTS ONE QUERY (A DNA SEQUENCE) AGAINST A DATABASE OF ONE SEQUENCE
        DB is the database to blast against
        (containing the sequence in the context of the plasmid)
        
        see BLAST parameters at
        https://biopython.org/docs/1.75/api/Bio.Blast.Applications.html
        
        -outfmt: 5 (xml), 10 (csv)
        --------------------------------------------------
        INPUTS:
        -- query [str]: a DNA sequence to test
        
        OUTPUTS:
        -- result [list]: a list of tuples
        [[(Acession,
           Evalue,
           query start,
           query end,
           target start,
           target end,
           (query strand, target strand))],
           [ ], ...]
        """
        
        db_file = '{}{}'.format(self.dir_blast_database, self.label)
        #database containing a single sequence
        query_file =  '{}BLAST_query_file_{}.txt'.format(self.dir_blast_tmp, self.label)
        xml_file = '{}BLAST_result_{}.xml'.format(self.dir_blast_tmp, self.label)
        
        #write query to a fasta file
        write_to_fasta_or_txt(query, query_file)
        
        #perform BLAST
        command = 'blastn -query {} -db {} -task blastn-short -perc_identity {} -evalue {} -word_size {} -reward {} -penalty {} -gapopen {} -gapextend {} -outfmt 5 -out {} -num_alignments 200 -num_threads 4'.format(query_file, db_file, self.PERC_IDENTITY, self.EVALUE, self.WORD_SIZE, self.REWARD, self.PENALTY, self.GAPOPEN, self.GAPEXTEND, xml_file)
        os.system(command)
        
        result_xml = open(xml_file)
        result_ncbixml = NCBIXML.read(result_xml)
        result_hsps = result_ncbixml.alignments[0].hsps
        result = [(result_hsps[i].score,
                   result_hsps[i].align_length,  #the length of the query sequence
                   result_hsps[i].query_start,
                   result_hsps[i].query_end,
                   result_hsps[i].sbjct_start,
                   result_hsps[i].sbjct_end,
                   result_hsps[i].strand,
                   result_hsps[i].query,
                   result_hsps[i].match,
                   result_hsps[i].sbjct,
                   result_hsps[i].positives/result_hsps[i].align_length) for i in range(len(result_hsps))]
       
        os.system('rm -f {}'.format(query_file)) #remove input file
        os.system('rm -f {}'.format(xml_file)) #remove output files
        
        #save BLAST result to both pickle and txt files
        with open('{}{}_BLAST_result.pickle'.format(self.dir_blast, self.label), 'wb') as f:
            pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)
        
        convert_blast_pickle_to_txt('{}{}_BLAST_result.pickle'.format(self.dir_blast, self.label),
                                    '{}{}_BLAST_result.txt'.format(self.dir_blast, self.label))
                                                 
        return result
    
    
    
    def calcu_blast_score(self, sequence):
        """
        Calculating BLAST matches to self, donor plasmids, and recipient plasmids
        """
        #sequence_database = [sequence] + [str(record.seq) for record in SeqIO.parse(self.plasmids_file, 'fasta')]
        #self.make_blast_database(sequence_database)
        
        #blast_result = self.make_blast_single(sequence)
                
        ##convert BLAST result to track (checked)
        #blast_result.pop(0) #remove self-identity (i.e., the first one). Then why you do blast for this one???
        #blast_score = [0] * len(sequence)
        #for k in range(len(blast_result)):
        #    pos_start = min(blast_result[k][2], blast_result[k][3]) - 1 #is it possible the other way? index started from 1
        #    pos_end = max(blast_result[k][2], blast_result[k][3]) - 1 #is it possible the other way? index started from 1
        #    blast_score[pos_start:pos_end+1] = [blast_result[k][0]] * (pos_end - pos_start + 1)
            
            
        sequence_database = [sequence] + [str(record.seq) for record in SeqIO.parse(self.plasmids_file, 'fasta')]
        self.make_blast_database(sequence_database)
        
        blast_result = self.make_blast_single(sequence)
        blast_result.pop(0) #remove self-identity (i.e., the first one).
        
        #convert BLAST result to track
        blast_score = [0] * len(sequence)
        if len(blast_result) > 0:
            for k in range(len(blast_result)):
                pos_start = blast_result[k][2] - 1 #index started from 1
                pos_end = blast_result[k][3] - 1 #index started from 1
                blast_score[pos_start:pos_end+1] = [blast_result[k][0]] * (pos_end - pos_start + 1)
            
        return blast_score
        
    
    
    def calcu_hairpin_score(self, sequence): #seems wrong tool for calculating hairpin!!!
        """
        """
        #window_size = self.hairpin_window_size
        #len_left = int(np.floor(window_size/2))
        #len_right = window_size - len_left - 1
        #sequence_pad = sequence[-len_left:] + sequence + sequence[:len_right]
        #sequence_slide_window = sliding_window(sequence_pad, window_size)
        #hairpin_score = [primer3.calc_hairpin(seq).dg for seq in sequence_slide_window]
        
        window_size = self.hairpin_window_size
        sequence_slide_window = sliding_window(sequence, window_size)
        tmp = int((window_size - 1)/2)
        hairpin_score = [0]*tmp + [primer3.calc_hairpin(seq).dg for seq in sequence_slide_window] + [0]*tmp
                
        return [-s for s in hairpin_score]
    
    
    
    #uncheck this
    def calcu_repeat_score(self, sequence):
        """
        """
        psAll = DU.getNDUAllSeq(sequence)
        idx_list = [k+1 for (k, ele) in enumerate(psAll)
                    if (ele >= self.NDU_threshold) and (k+1) < self.NDU_repeat_window_size]
        #limit to repeats that are above threshold and that are spaced under repeat window size

        repeat_score = [0] * len(sequence)
        for k in idx_list:
            tmp_1 = DU.getNDUWindowSeq(sequence, k, self.NDU_repeat_window_size)
            tmp_2 = [0] * int(self.NDU_repeat_window_size/2)
            tmp_3 = tmp_2 + tmp_1 + tmp_2
            tmp_4 = [x if x >= self.NDU_threshold else 0 for x in tmp_3]
            #remove noise, consider only values above threshold
            repeat_score = [x1-x2 for (x1, x2) in zip(repeat_score, tmp_4)]
            #position score vector so far
  
        return [-s for s in repeat_score]
        
        
    #unfinished
    #def calcu_mfe_score(self, sequence):
    #    return
    
        
    def gc_conent_to_score(self, x):
        """
        CONVERT GC CONTENT TO GC SCORE
        -----------------------------------------
        GC score is a value in the range of [0, 100],
        with 0 being the best, and 100 being the worst.
        """
        GC_threshold_lower_limit = 0.2
        GC_threshold_upper_limit = 0.8
        score_scale = 10
        
        if x <= GC_threshold_lower_limit:
            return score_scale
        
        elif x <= self.GC_threshold_lower:
            y = (x - self.GC_threshold_lower)**2 / (GC_threshold_lower_limit - self.GC_threshold_lower)**2 * score_scale
            return y
            
        elif x <= self.GC_threshold_upper:
            return 0
        
        elif x <= GC_threshold_upper_limit:
            y = (x - self.GC_threshold_upper)**2 / (GC_threshold_upper_limit - self.GC_threshold_upper)**2 * score_scale
            return y
        else:
            return score_scale
        
        
    
    def norm_score(self, score, option):
        """
        """
        if option == 'gc':
            #change GC content to a normalized GC score:
            #score_normalized = [i*10 if i < self.GC_threshold_lower or i > self.GC_threshold_upper else 0 for i in score]
            score_normalized = [self.gc_conent_to_score(x) for x in score]
        
        elif option == 'blast':
            #score_normalized = [s/2 for s in score]
            score_normalized = [s/10 for s in score]
        
        elif option == 'hairpin':
            score_normalized = [s/500 if s>self.hairpin_threshold else 0 for s in score]
        
        elif option == 'repeat':
            score_normalized = score
            
        #elif option == 'mfe':
                
        return score_normalized
    

    

    def calcu_score(self, sequence):
        """
        """
        score_dict = dict()
        if self.use_GC:
            my_print('      1) GC score...', not self.parallel)
            gc_content = self.calcu_gc_content(sequence)
            
            gc_score_normalized = self.norm_score(gc_content, 'gc')
            score_dict['gc_score'] = {'raw': gc_content,
                                      'norm': gc_score_normalized}
            
        if self.use_BLAST:
            my_print('      2) BLAST score...', not self.parallel)
            blast_score = self.calcu_blast_score(sequence)
            blast_score_normalized = self.norm_score(blast_score, 'blast')
            score_dict['blast_score'] = {'raw': blast_score,
                                         'norm': blast_score_normalized}
            
        if self.use_Primer3:
            my_print('      3) Primer3 hairpin score...', not self.parallel)
            hairpin_score = self.calcu_hairpin_score(sequence)
            hairpin_score_normalized = self.norm_score(hairpin_score, 'hairpin')
            score_dict['primer3_hairpin_score'] = {'raw': hairpin_score,
                                                   'norm': hairpin_score_normalized}
              
            #print('aaa---')
            #print(hairpin_score)
            
        if self.use_NDU:
            my_print('      4) NDU repeat score...', not self.parallel)
            repeat_score = self.calcu_repeat_score(sequence)
            repeat_score_normalized = self.norm_score(repeat_score, 'repeat')
            score_dict['ndu_repeat_score'] = {'raw': repeat_score,
                                              'norm': repeat_score}
            
        
        #if self.use_MFE:

            
        score_normalized = np.array([score_dict[ele]['norm']
                                     for ele in score_dict.keys()],
                                    dtype=float)
        score_dict['summed_score'] = list(np.sum(score_normalized, axis=0))
        
        return score_dict
        
        
        
def score(sequence, param_dict, label):
    """
    """
    my_object = CalculateScore(sequence, param_dict, label)
    result = my_object.calcu_score(sequence)
    
    return result



