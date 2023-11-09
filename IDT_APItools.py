#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 11:22:14 2022
Functions for testing if sequences can be synthesized at IDT
Uses the IDT API
To use the API, an access token must be generated first using get_access_token
This token lasts ~1h and can be used for multiple queries
Other functions are specific for different IDT processes (Gblocks, Eblocks, etc.)
@author: sashalevy
"""

################## FUNCTIONS #######################

from __future__ import print_function
from base64 import b64encode
import json
from urllib import request, parse
import http.client


def get_access_token(client_id, client_secret, idt_username, idt_password):
    """
    Create the HTTP request, transmit it, and then parse the response for the 
    access token.
    
    The body_dict will also contain the fields "expires_in" that provides the 
    time window the token is valid for (in seconds) and "token_type".
    """

    # Construct the HTTP request
    authorization_string = b64encode(bytes(client_id + ":" + client_secret, "utf-8")).decode()
    request_headers = { "Content-Type" : "application/x-www-form-urlencoded",
                        "Authorization" : "Basic " + authorization_string }
                    
    data_dict = {   "grant_type" : "password",
                    "scope" : "test",
                    "username" : idt_username,
                    "password" : idt_password }
    request_data = parse.urlencode(data_dict).encode()

    post_request = request.Request("https://www.idtdna.com/Identityserver/connect/token", 
                                    data = request_data, 
                                    headers = request_headers,
                                    method = "POST")

    # Transmit the HTTP request and get HTTP response
    response = request.urlopen(post_request)

    # Process the HTTP response for the desired data
    body = response.read().decode()
    
    # Error and return the response from the endpoint if there was a problem
    if (response.status != 200):
        raise RuntimeError("Request failed with error code:" + response.status + "\nBody:\n" + body)
    
    body_dict = json.loads(body)
    return body_dict["access_token"]




def Gblock_complexity(seq, access_token):
    """
    Determine the Gblock complexity score and areas that may need to be recoded
    Scores <7 are acceptable
    Outputs a dictrionary with the data
    Output keys: 
        'ActualValue' 
        'DisplayText'
        'ForwardLocations'
        'IsViolated' - T/F on whether there is a potentially problematic regions
        'Name'
        'RepeatedSegment'
        'ReverseLocations'
        'Score' - the complexity score (lower is better)
        'StartIndex'
        'TerminalEnd'
        'ThresholdOutput': e.g. {'Value': 40, 'ThresholdType': 1, 'MinLength': 0, 'MaxLength': 0, 'MinPercentage': 0.0, 'MaxPercentage': 0.0, 'WindowLength': 0, 'Quantity': 0} 
        'ServiceProductId'
        'MinimumRepeatLength'
        'RepeatPercentage'
        'GCPercentage'
        'Length'
        'Rank'
    
    """
    conn = http.client.HTTPSConnection("www.idtdna.com")
    payload = "[\n    {\n        \"Name\": \"My gBlock2\",\n        \"Sequence\":\"" + seq + "\n\n\"}]"
    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer ' + access_token 
        }
    conn.request("POST", "/api/complexities/screengBlockSequences", payload, headers)
    res = conn.getresponse()
    data = res.read()
    pdata = data.decode("utf-8") #python readible data
    if pdata != '[[]]': #has some poblematic DNA sequence
        pdata = pdata[1:-1] #remove first and last 3 characters
        ddata = json.loads(pdata) # a dictionary with the output data
    else: #has no problematic DNA sequence
        ddata = {'ActualValue': 0.0, 'DisplayText': 'None', 'ForwardLocations': [], 'IsViolated': False, 'Name': 'Overall Repeat', 'RepeatedSegment': None, 'ReverseLocations': [], 'Score': 0.0, 'StartIndex': 0, 'TerminalEnd': 0, 'ThresholdOutput': {'Value': 0, 'ThresholdType': 1, 'MinLength': 0, 'MaxLength': 0, 'MinPercentage': 0.0, 'MaxPercentage': 0.0, 'WindowLength': 0, 'Quantity': 0}, 'ServiceProductId': 0, 'MinimumRepeatLength': 0, 'RepeatPercentage': 0.0, 'GCPercentage': 0.0, 'Length': 0, 'Rank': 0.0}
    return(ddata)


def Gblock_hifi_complexity(seq, access_token):
    """
    Determine the GblockHiFi complexity score and areas that may need to be recoded
    Scores <7 are acceptable
    Outputs a dictrionary with the data
    Output keys: 
        'ActualValue' 
        'DisplayText'
        'ForwardLocations'
        'IsViolated' - T/F on whether there is a potentially problematic regions
        'Name'
        'RepeatedSegment'
        'ReverseLocations'
        'Score' - the complexity score (lower is better)
        'StartIndex'
        'TerminalEnd'
        'ThresholdOutput': {'Value': 40, 'ThresholdType': 1, 'MinLength': 0, 'MaxLength': 0, 'MinPercentage': 0.0, 'MaxPercentage': 0.0, 'WindowLength': 0, 'Quantity': 0} 
        'ServiceProductId'
        'MinimumRepeatLength'
        'RepeatPercentage'
        'GCPercentage'
        'Length'
        'Rank'
    """
    conn = http.client.HTTPSConnection("www.idtdna.com")
    payload = "[\n    {\n        \"Name\": \"My gBlock2\",\n        \"Sequence\":\"" + seq + "\n\n\"}]"
    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer ' + access_token 
        }
    conn.request("POST", "/restapi/v1/Complexities/ScreenGBlockHifiSequences", payload, headers)
    res = conn.getresponse()
    data = res.read()
    pdata = data.decode("utf-8") #python readible data
    if pdata != '[[]]': #has some poblematic DNA sequence
        pdata = pdata[1:-1] #remove first and last 3 characters
        ddata = json.loads(pdata) # a dictionary with the output data
    else: #has no problematic DNA sequence
        ddata = {'ActualValue': 0.0, 'DisplayText': 'None', 'ForwardLocations': [], 'IsViolated': False, 'Name': 'Overall Repeat', 'RepeatedSegment': None, 'ReverseLocations': [], 'Score': 0.0, 'StartIndex': 0, 'TerminalEnd': 0, 'ThresholdOutput': {'Value': 0, 'ThresholdType': 1, 'MinLength': 0, 'MaxLength': 0, 'MinPercentage': 0.0, 'MaxPercentage': 0.0, 'WindowLength': 0, 'Quantity': 0}, 'ServiceProductId': 0, 'MinimumRepeatLength': 0, 'RepeatPercentage': 0.0, 'GCPercentage': 0.0, 'Length': 0, 'Rank': 0.0}
    return(ddata)

def Eblock_complexity(seq, access_token):
    """
    Determine the Eblock complexity score and areas that may need to be recoded
    Scores <15 are acceptable
    Max length = 900
    Null data means no significant complexity
    """
    conn = http.client.HTTPSConnection("www.idtdna.com")
    payload = "[\n    {\n        \"Name\": \"My gBlock2\",\n        \"Sequence\":\"" + seq + "\n\n\"}]"
    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer ' + access_token 
        }
    conn.request("POST", "/restapi/v1/Complexities/ScreenEBlockSequences", payload, headers)
    res = conn.getresponse()
    data = res.read()
    pdata = data.decode("utf-8") #python readible data
    if pdata != '[[]]': #has some poblematic DNA sequence
        pdata = pdata[1:-1] #remove first and last 3 characters
        ddata = json.loads(pdata) # a dictionary with the output data
    else: #has no problematic DNA sequence
        ddata = {'ActualValue': 0.0, 'DisplayText': 'None', 'ForwardLocations': [], 'IsViolated': False, 'Name': 'Overall Repeat', 'RepeatedSegment': None, 'ReverseLocations': [], 'Score': 0.0, 'StartIndex': 0, 'TerminalEnd': 0, 'ThresholdOutput': {'Value': 0, 'ThresholdType': 1, 'MinLength': 0, 'MaxLength': 0, 'MinPercentage': 0.0, 'MaxPercentage': 0.0, 'WindowLength': 0, 'Quantity': 0}, 'ServiceProductId': 0, 'MinimumRepeatLength': 0, 'RepeatPercentage': 0.0, 'GCPercentage': 0.0, 'Length': 0, 'Rank': 0.0}
    return(ddata)

def Gene_complexity(seq, access_token):
    """
    Determine the Gene complexity score and areas that may need to be recoded
    Scores <7 are acceptable
    Max length = 5000
    Outputs a dictrionary with the data
    Output keys: 
        'ActualValue' 
        'DisplayText'
        'ForwardLocations'
        'IsViolated' - T/F on whether there is a potentially problematic regions
        'Name'
        'RepeatedSegment'
        'ReverseLocations'
        'Score' - the complexity score (lower is better)
        'StartIndex'
        'TerminalEnd'
        'ThresholdOutput': {'Value': 40, 'ThresholdType': 1, 'MinLength': 0, 'MaxLength': 0, 'MinPercentage': 0.0, 'MaxPercentage': 0.0, 'WindowLength': 0, 'Quantity': 0} 
        'ServiceProductId'
        'MinimumRepeatLength'
        'RepeatPercentage'
        'GCPercentage'
        'Length'
        'Rank'
    """
    conn = http.client.HTTPSConnection("www.idtdna.com")
    payload = "[\n    {\n        \"Name\": \"My gBlock2\",\n        \"Sequence\":\"" + seq + "\n\n\"}]"
    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer ' + access_token 
        }
    conn.request("POST", "/restapi/v1/Complexities/ScreenGeneSequences", payload, headers)
    res = conn.getresponse()
    data = res.read()
    pdata = data.decode("utf-8") #python readible data
    if pdata != '[[]]': #has some poblematic DNA sequence
        pdata = pdata[1:-1] #remove first and last 3 characters
        ddata = json.loads(pdata) # a dictionary with the output data
    else: #has no problematic DNA sequence
        ddata = {'ActualValue': 0.0, 'DisplayText': 'None', 'ForwardLocations': [], 'IsViolated': False, 'Name': 'Overall Repeat', 'RepeatedSegment': None, 'ReverseLocations': [], 'Score': 0.0, 'StartIndex': 0, 'TerminalEnd': 0, 'ThresholdOutput': {'Value': 0, 'ThresholdType': 1, 'MinLength': 0, 'MaxLength': 0, 'MinPercentage': 0.0, 'MaxPercentage': 0.0, 'WindowLength': 0, 'Quantity': 0}, 'ServiceProductId': 0, 'MinimumRepeatLength': 0, 'RepeatPercentage': 0.0, 'GCPercentage': 0.0, 'Length': 0, 'Rank': 0.0}
    return(ddata)


