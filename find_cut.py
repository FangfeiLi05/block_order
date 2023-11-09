import pandas as pd
import numpy as np
import time
import argparse
from sub_functions import my_print

class CutBlocks(object):
    def __init__(self, position_score, n, l_min, l_max, o_min, o_max):
        self.position_score = position_score
        self.n = n #number of fragments
        self.l_min, self.l_max = l_min, l_max
        self.o_min, self.o_max = o_min, o_max
        self.L = len(position_score)
        
        #if self.l_min <= self.o_max:
        #    print('break beacuse: l_min > o_max not satisfied')
        #if self.n_min > self.n_max:
        #    print('break because: no valid n found')
        
        self.opt_obj = 'max'
       

    def startState(self):
        return (1, 1, self.n)
        
        
    def isEnd(self, state_s, state_e, state_k):
        return (state_s, state_e, state_k) == (self.L, self.L, 0)
        
        
    def calcuCost(self, state_s, state_e, state_k):
        # return cost of a state
        if not self.isEnd(state_s, state_e, state_k):
            return sum(self.position_score[state_s-1: state_e])
        else:
            return 0

    
    def succAndCost(self, state_s, state_e, state_k):
        # return a list of (state_s, state_e, state_k, cost)
        
        if state_k != 1:
            skip_size_list = [20, 10, 5, 2, 1]
            #skip_size_list = [10, 5, 2, 1]
            #skip_size_list = [1]
            result = []
            i = 0
            while len(result) == 0:
                new_e_list = [new_e for new_e in range(state_s + self.l_min - 1,
                                                       state_s + self.l_max - 1 + 1,
                                                       skip_size_list[i])]

                b1_b2_list = [(max(new_e - self.o_max,
                                   self.L - (state_k - 1) * self.l_max + (state_k - 2) * self.o_min) + 1,
                               min(new_e - self.o_min,
                                   self.L - (state_k - 1) * self.l_min + (state_k - 2) * self.o_max) + 1)
                               for new_e in new_e_list]
            
                result = [(new_s, new_e, state_k-1, self.calcuCost(new_s, new_e, state_k-1))
                          for (j, new_e) in enumerate(new_e_list)
                          if b1_b2_list[j][1] - b1_b2_list[j][0] >= 0
                          for new_s in range(b1_b2_list[j][0], b1_b2_list[j][1] + 1)]
                i += 1
                
            return result
        
        else:
            return [(self.L, self.L, 0, 0)]
            

### Algorithms
def dynamicProgramming(problem):
    cache = {} # state -> futureCost(state)
    def futureCost(state_s, state_e, state_k):
        # Base case
        if problem.isEnd(state_s, state_e, state_k):
            return 0
            
        tmp = '{}_{}_{}'.format(state_s, state_e, state_k)
        if tmp in cache: # Exponential savings
            return cache[tmp][0]
            
        # Actually doing work
        if problem.opt_obj == 'sum':
            result = min((cost+futureCost(newState_s, newState_e, newState_k),
                          newState_s,
                          newState_e,
                          newState_k,
                          cost)
                         for newState_s, newState_e, newState_k, cost
                         in problem.succAndCost(state_s, state_e, state_k))
            
        if problem.opt_obj == 'max':
            result = min((max(cost, futureCost(newState_s, newState_e, newState_k)),
                          newState_s,
                          newState_e,
                          newState_k,
                          cost)
                         for newState_s, newState_e, newState_k, cost
                         in problem.succAndCost(state_s, state_e, state_k))

        cache[tmp] = result
        return result[0]
        
    (state_s, state_e, state_k) = problem.startState()
    totalCost = futureCost(state_s, state_e, state_k)
       
   
    # Recover history
    history = []
    while not problem.isEnd(state_s, state_e, state_k):
        tmp = '{}_{}_{}'.format(state_s, state_e, state_k)
        _, newState_s, newState_e, newState_k, cost = cache[tmp]
        history.append((newState_s, newState_e, newState_k, cost))
        state_s, state_e, state_k = (newState_s, newState_e, newState_k)
    
    return (totalCost, history)
    


def cut_fixed_n(position_score, n, l_min=1000, l_max=1744, o_min=50, o_max=50):
    """
    """
    L = len(position_score)
    L_part = 100000
    
    if L <= L_part:
        problem = CutBlocks(position_score, n, l_min, l_max, o_min, o_max)
        totalCost, history = dynamicProgramming(problem)

        x_frag_list = [1] + [history[i][0] for i in range(len(history)-1)]
        y_frag_list = [history[i][1] for i in range(len(history))]
        
    elif L <= 2*L_part: # not finish yet
        problem = CutBlocks(position_score[:L_part], n)
        print('\nRunning 1st DNA cuts algorithm ...')
        totalCost, history = dynamicProgramming(problem)

        x_frag_list = [1] + [history[i][0] for i in range(len(history)-1)]
        y_frag_list = [history[i][1] for i in range(len(history))]
        len_frag_list = [y_frag_list[i] - x_frag_list[i] for i in range(len(x_frag_list))]
        
        tmp = [(k, ele) for (k, ele) in enumerate(x_frag_list) if (L - ele) < L_part][0] ###
        problem = CutBlocks(position_score[tmp[1]-1:])
        print('\nRunning 2nd DNA cuts algorithm ...')
        totalCost, history = dynamicProgramming(problem)
        
        x_frag_list = x_frag_list[:tmp[0]] + [1 + tmp[1] - 1] + [history[i][0] + tmp[1] for i in range(len(history)-1)]
        y_frag_list = y_frag_list[:tmp[0]] + [history[i][1] + tmp[1] - 1 for i in range(len(history))]
    
    return (totalCost, x_frag_list, y_frag_list)
    
    
def cut(position_score, n_max=10, l_min=1000, l_max=1744, o_min=50, o_max=50, print_option=True):
    """
    """
    t_start = time.time()
        
    L = len(position_score)
    n_min = int(np.ceil((L - o_min) / (l_max - o_min)))
    n_max = np.min([int(np.floor((L - o_max) / (l_min - o_max))), n_max])
    
    my_print('      1) parameters:'.format(L), print_option)
    my_print('         total length: {}'.format(L), print_option)
    my_print('         length for a fragment: [{}, {}]'.format(l_min, l_max), print_option)
    my_print('         overlap between two fragments: [{}, {}]'.format(o_min, o_max), print_option)
        
    output_dict = dict()
    for n in range(n_min, n_max+1):
        output_dict[n] = cut_fixed_n(position_score, n, l_min, l_max, o_min, o_max)
                
    (totalCost_opt, n_opt, x_frag_opt, y_frag_opt) = min((output_dict[i][0],
                                                          i,
                                                          output_dict[i][1],
                                                          output_dict[i][2]) for i in output_dict.keys())
    
    len_frag_opt = [y_frag_opt[i] - x_frag_opt[i] for i in range(len(x_frag_opt))]
         
    t_end = time.time()
    
    my_print('      2) optimized cost: {}'.format(totalCost_opt), print_option)
    my_print('      3) result:', print_option)
    for i in range(len(x_frag_opt)):
        my_print('         {}: pos=[{}, {}], len={}'.format(i+1, x_frag_opt[i], y_frag_opt[i], len_frag_opt[i]), print_option)
        

    result_dict = {'optimal_cost': totalCost_opt,
                   'block_number': n_opt,
                   'start_pos': x_frag_opt,
                   'end_pos': y_frag_opt}
    
    return result_dict
    
    
    


def main():
    """
    """
    parser = argparse.ArgumentParser(description='Find optimal DNA fragments',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-s', '--score', type=str, required=True,
                        help='a .csv file with a column being the position score for each base pair')
    # The position score for each bp is non-negative. The smaller the score, the better.
                    
    parser.add_argument('-fl', '--fragment_length', nargs='*', type=int, default = (1000,1744),
                        help='minimum and maximum length allowed for a fragment')
    
    parser.add_argument('-ol', '--overlap_length', nargs='*', type=int, default = (50,50),
                        help='minimum and maximum overlapping length allowed between two fragments')
                        
    parser.add_argument('-nm', '--number_max', nargs='*', type=int, default = 10,
                        help='maximum number of fragments allowed')
        
    args = parser.parse_args()
    
    t_start = time.time()
        
    #position_score = list(np.array(pd.read_csv(args.score, header=None)[0], dtype=float)*(-1))
    position_score = list(pd.read_csv(args.score, header=None)[0])
    L = len(position_score)
    
    l_min, l_max = args.fragment_length
    o_min, o_max = args.overlap_length
    n_max = args.number_max
    
    n_min = int(np.ceil((L - o_min) / (l_max - o_min)))
    #n_max = int(np.floor((L - o_max) / (l_min - o_max)))
    n_max = np.min([int(np.floor((L - o_max) / (l_min - o_max))), n_max])
    
    print('\nRunning DNA fragments cutting algorithm ...')
    print('length: {}'.format(L))
    print('l_min: {}, l_max: {}'.format(l_min, l_max))
    print('o_min: {}, o_max: {}'.format(o_min, o_max))
    print('n_min: {}, n_max: {}'.format(n_min, n_max))
    
    output_dict = dict()
    for n in range(n_min, n_max+1):
        output_dict[n] = cut_fixed_n(position_score, n, l_min, l_max, o_min, o_max)
                
    (totalCost_opt, n_opt, x_frag_opt, y_frag_opt) = min((output_dict[i][0], i, output_dict[i][1], output_dict[i][2])
                                                         for i in output_dict.keys())
    
    len_frag_opt = [y_frag_opt[i] - x_frag_opt[i] for i in range(len(x_frag_opt))]
            
    print('-- total cost: {}'.format(totalCost_opt))
    print('-- number of fragments: {}'.format(n_opt))
    print('-- for each fragment ...')
    for i in range(len(x_frag_opt)):
        print('---- frag_{}: {} ... {}, len={}'.format(i+1, x_frag_opt[i], y_frag_opt[i], len_frag_opt[i]))
    
    t_end = time.time()
    print('-- running time: {}'.format(t_end - t_start))
    print('Done!')
    
    
if __name__ == "__main__":
    main()

    
