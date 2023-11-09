#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 12:20:13 2022
Functions for testing if sequences can be synthesized at Twist
Uses the Twist API
@author: sashalevy

Developer's Portal:
https://developers.twistdna.com/docs/tapi

API CREDENTIALS

Authorization:
JWT ***

X-End-User-Token: ***

Sandbox (Staging):https://twist-api.twistbioscience-staging.com
Production:https://twist-api.twistdna.com
"""

def submit_gene_block(seq, name = "TEST_FRAGMENT", adaptors = False ):
    """
    Create a new construct
    Construct creation for Genes and Oligo Pools invokes an asynchronous scoring process.
    There are 2 options in order to receive the results of the scoring process:
        Pass a callback_url query param (e.g POST /v1/users/{email}/constructs/?callback_url=...) - Scoring results will be sent to the provided URL, possibly in multiple requests.
        Poll for results at the Scoring Endpoint
        sequences - Must be an array for OLIGO_POOL, but for other types an object is required for each sequence.
        
        Length of each sequence must fit the limitations of the construct type and must contain only upper case ACTG.
        
        You can find the information about the vectors available for your use, at GET /v1/users/{email}/vectors/ endpoint.
        
        If the construct type is "CLONED_GENE": - insertion_point_mes_uid - must be a valid insertion site ID. Make sure insertion_point_mes_uid suits the length of your sequence.
        If your account is enabled for extra long clonal genes (up to 5000bp):
            Mixed orders of genes longer than 3350bp and shorter than 3350bp are not allowed.
            Order size for genes longer than 3350bp must be limited to 500.
            IMPORTANT: Throttling

    We are currently limiting each end user to a specific number of simultaneous scoring requests.

    In case you get an error, you can contact us directly so we can increase your limits if you require it.
    
    Returns a dictionary with the construct id and other info
    
    """
    import requests
    import json
    url = "***"
    payload = {
        "sequences": seq,
        "name": name,
        "type": "NON_CLONED_GENE",
        "adapters_on": adaptors
    }
    headers = {
        "Content-Type": "application/json",
        "X-End-User-Token": "***",
        "Authorization": "JWT *** "
    }
    response = requests.request("POST", url, json=payload, headers=headers)
    data = json.loads(response.text) # a dictionary with the output data
    return(data)
    


    
def get_construct_score(construct):
    """
    Use this API to poll a list of constructs objects and their Scoring results.
    Sequences will not be sent back and cannot be modified.
    Unscored constructs can not be edited or used in a quote/order
    Errors:
    - 4000 - General error
    - 4001 - Problematic sequence
    - 4002 - Repeats or extreme high/low GC in the highlighted region may have caused this problem
    - 4003 - Provided sequence is invalid
    - 4004 - Provided sequence contains invalid character(s)
    - 4100 - Invalid sequence length
    - 4101 - Sequence is too short
    - 4102 - Sequence is too long
    - 4103 - Sequences longer than 1,700 bp elevate risk marginally. Consider splitting your sequence into smaller pieces to increase the likelihood of success
    - 4104 - Sequences longer than 3,100 bp elevate risk. Consider splitting your sequence into smaller pieces to increase the likelihood of success
    - 4200 - Invalid GC content
    - 4201 - The overall GC content of your sequence must be under 65%. Under 60% will be optimal for success
    - 4202 - The overall GC content of your sequence must be over 25%
    - 4203 - The difference in GC content between the highest-GC and lowest-GC 50bp windows exceeds 52%. Please even out the peaks/troughs in GC content of your sequence. Lowering the peaks will be optimal for success
    - 4300 - Secondary structure
    - 4301 - Hairpin detected
    - 4303 - Long direct repeat (greater or equal to 20bp) detected. Please break up repeats. Fewer/lower homology repeats will be optimal for success
    - 4304 - Direct Repeat with high Tm (greather than 60C) detected. Please break up repeats. Fewer/lower homology repeats will be optimal for success
    - 4305 - More than 45% of your sequence is composed of small repeats (9bp or longer). This increases complexity. Please break up repeats, perhaps by varying your codon usage
    - 4306 - Long, low homology repeat region detected. Please break up repeats. Fewer/lower homology repeats will be optimal for success.
    - 4401 - Attempted hierarchical design but could not find an acceptable split point
    - 4402 - We are unable to make the sequence as is. Please try to run our Codon Optimization or try to split the gene into two parts
    - 4403 - Attempted hierarchical design but could not find an acceptable split point
    - 4404 - Design request results in fragment sizes less than 300bp; please retry design with fewer fragments
    - 4405 - Design request results in fragment sizes greater than 1,800bp; please retry design with more fragments
    - 4501 - His tags consisting of 5 or more identical codons increase complexity. Please vary the codons used (e.g. CAT CAC CAT... instead of CAT CAT CAT...)
    - 4502 - CpG multimeric segments of 14 or more bases increase complexity. Please break up these low-complexity sequences
    - 4503 - Long homopolymer stretches increase complexity. Please break up long homopolymers
    - 4504 - Sequence contains one or more Gateway cloning att sites
    - 4505 - Sequence contains one or more impermissible sequences that complicate manufacture or QC.
    - 4506 - Clonal gene contains sub-sequence with high homology to CCDB.
    - 4507 - We detected a restriction site of a methylation sensitive enzyme that overlaps with a known methylation site.
    - 4508 - Unable to design primers for this construct. Please increase the GC content in the first and last 60 bases of your sequence.
    - 5001 - Unable to split fragment for synthesis.
    - 5002 - Fragment design exception
    """
    import requests
    import json
    construct_id = construct["id"]
    url = "***"
    querystring = {"id__in": construct_id, "scored":"true" }
    headers = {
        "Content-Type": "application/json",
        "X-End-User-Token": "***",
        "Authorization": "JWT ***"
    }
    data = {'key':'value'}
    while type(data) == dict: #changes to a list when the score has been generated
        response = requests.request("GET", url, headers=headers, params=querystring)
        data = json.loads(response.text) # a dictionary with the output data    
    return(data[0])


def submit_codon_optimization(construct): 
    '''
    Get previously created Codon Fitting Optimization Jobs based on their construct IDs.
    To create a new Codon Optimization job, please use the POST Codon Optimization endpoint
    '''
    
    import requests
    url = "https://twist-api.twistbioscience-staging.com/v1/users/sflevy@stanford.edu/codon-optimizations/"
    payload = [
        {
            "sequence": seq,
            "external_id": "test_codon_opt",
            "organism": "Escherichia coli",
            "avoid_introducing": ["EcoRI", "XbaI"],
            "optimization_start": 3,
            "optimization_len": 9,
            "adapters_on": True,
            "old_scoring": False
        }
    ]
    headers = {
        "Content-Type": "application/json",
        "X-End-User-Token": "***",
        "Authorization": "JWT ***"
    }
    response = requests.request("POST", url, json=payload, headers=headers)
    print(response.text)

#submit_codon_optimization(natR)

def get_codon_optimization(seq): 
    import requests
    url = "***"
    querystring = {"id__in":"***","completed":"true"}
    headers = {
        "Content-Type": "application/json",
        "X-End-User-Token": "***",
        "Authorization": "JWT ***"
    }
    response = requests.request("GET", url, headers=headers, params=querystring)
    print(response.text)

