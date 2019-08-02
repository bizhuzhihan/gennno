#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Gene annotation
# Copyright (c) 2019 Zhihan Zhu 
# https://github.com/bizhuzhihan/gennno

import re
import time
import subprocess
import pandas as pd
from urllib import request
from bs4 import BeautifulSoup as bs

welcome = b'\n\xe6\x94\xaf\xe6\x8c\x81\xe7\x9a\x84: CARD, UniProt\n\xe4\xbd\xbf\xe7\x94\xa8\xe8\xaf\xb4\xe6\x98\x8e: gennno.help()\n\n\xe2\x88\xa0( \xe1\x90\x9b \xe3\x80\x8d\xe2\x88\xa0)\xef\xbc\xbf \xe2\x88\xa0( \xe1\x90\x9b \xe3\x80\x8d\xe2\x88\xa0)\xef\xbc\xbf \xe2\x88\xa0( \xe1\x90\x9b \xe3\x80\x8d\xe2\x88\xa0)\xef\xbc\xbf \xe2\x88\xa0( \xe1\x90\x9b \xe3\x80\x8d\xe2\x88\xa0)\xef\xbc\xbf \xe2\x88\xa0( \xe1\x90\x9b \xe3\x80\x8d\xe2\x88\xa0)\xef\xbc\xbf \xe2\x88\xa0( \xe1\x90\x9b \xe3\x80\x8d\xe2\x88\xa0)\xef\xbc\xbf  \n    '
print(welcome.decode('utf-8'))

####################################################################################################
# Class
class dove:
    def __init__(self, genes, path):
        self.genes = genes
        self.db = pd.DataFrame()
        self.db['gene'] = self.genes
        self.path = path

        self.b_card = 0
        self.b_uniprot = 0
        self.b_string = 0

    ################################################################################################
    # CARD Database
    def card(self, path_CARD=None):
        if(self.b_card == 0):
            if(path_CARD==None):
                request_url = "https://card.mcmaster.ca/latest/data"
                request.urlretrieve(request_url, self.path+"card.tar.bz2")
                subprocess.run(['tar', '-xjvf', self.path+"card.tar.bz2", '-C', self.path])
                df = pd.read_json(self.path+'card.json').T.sort_values('ARO_name')
                df['gene'] = df['ARO_name']
            else:
                df = pd.read_json(path_CARD).T.sort_values('ARO_name')
                df['gene'] = df['ARO_name']
        
            self.db = pd.merge(self.db, df[['gene','ARO_description']], on='gene', how='left')
            self.b_card = 1

        
    ################################################################################################
    # Uniprot Database
    def uniprot(self, prior):
        if(self.b_uniprot ==0):
            tab_head = "https://www.uniprot.org/uniprot/?query=gene_exact:"
            tab_tail = "&format=tab&force=true&columns=id,protein%20names,organism&sort=score&compress=no"
            uniprot_db = "https://www.uniprot.org/uniprot/"
            uniprot_gene = []
            uniprot_function = []

            for i in range(len(self.genes)):
                request_url = tab_head + self.genes[i] + tab_tail
                df = pd.read_table(request_url)
                time.sleep(1)
                for j in range(df.shape[0]):
                    try:
                        # a=((re.search(prior, df['Organism'][j], flags=0))!=None)
                        # b=((re.search(product[i], df['Protein names'][j], flags=0))!=None)
                        # if(a&b):
                        if(prior==df['Organism'][j]):
                            with request.urlopen(uniprot_db+df['Entry'][j]+'.xml') as f:
                                data = f.read()
                                data = str(data)
                            time.sleep(1)
                            bb = bs(data, "lxml")
                            bb = bb.find('comment', type="function")
                            try:
                                bb = bb.get_text()
                                bb = re.sub(r'\\n', "", bb)
                                uniprot_gene.append(self.genes[i])
                                uniprot_function.append(bb)
                                break
                            except:
                                break
                    except:
                        continue

            df = pd.DataFrame()
            df['gene'] = uniprot_gene
            df['uniprot_function'] = uniprot_function
            self.db = pd.merge(self.db, df, on='gene', how='left')
            self.b_uniprot = 1

    
    # ################################################################################################
    # # STRING Database
    # def string(self, species, limit):
    #     if(self.b_string == 0):
    #         string_api_url = "https://string-db.org/api"
    #         output_format = "json"
    #         method = "interaction_partners"

    #         request_url = string_api_url + "/" + output_format + "/" + method + "?"
    #         request_url += "identifiers=%s" % "%0d".join(self.genes)
    #         request_url += "&" + "species=" + species
    #         request_url += "&" + "limit=" + str(limit)
    #         # data = pd.read_json(request_url)
    #         # hhhhhh
    #         # 
    #         # self.db = pd.merge(self.db, df, on='gene', how='left')
    #         # self.b_string = 1

####################################################################################################
# Help
def help():
    doc = b'\nCARD:\neg = gennno.dove(genelist, path)\neg.card(path_card)\neg.db\n# \xe5\x8f\x82\xe6\x95\xb0\xe8\xaf\xb4\xe6\x98\x8e\ngenelist: \xe5\x9f\xba\xe5\x9b\xa0\xe5\x90\x8d\xe7\xa7\xb0\xef\xbc\x88list\xef\xbc\x89\npath: \xe5\xb7\xa5\xe4\xbd\x9c\xe8\xb7\xaf\xe5\xbe\x84\xef\xbc\x8c\xe7\x94\xa8\xe4\xba\x8e\xe4\xbf\x9d\xe5\xad\x98\xe4\xb8\xb4\xe6\x97\xb6\xe6\x96\x87\xe4\xbb\xb6\npath_card: card.json\xe7\x9a\x84\xe8\xb7\xaf\xe5\xbe\x84\xef\xbc\x8c\xe6\xb2\xa1\xe6\x9c\x89\xe5\xb0\xb1\xe4\xb8\x8d\xe5\xa1\xab\xef\xbc\x8c\xe4\xbc\x9a\xe8\x87\xaa\xe5\x8a\xa8\xe4\xb8\x8b\xe8\xbd\xbd\xe5\x88\xb0\xe5\xb7\xa5\xe4\xbd\x9c\xe8\xb7\xaf\xe5\xbe\x84\xe4\xb8\xad\neg.db: \xe8\xbf\x94\xe5\x9b\x9e\xe5\x8c\x85\xe5\x90\xabgene\xef\xbc\x8cdescription\xe7\x9a\x84pd.DataFrame\n\nUniprot:\neg = gennno.dove(genelist, path)\neg.uniprot(species_name)\neg.db\n# \xe5\x8f\x82\xe6\x95\xb0\xe8\xaf\xb4\xe6\x98\x8e\ngenelist: \xe5\x9f\xba\xe5\x9b\xa0\xe5\x90\x8d\xe7\xa7\xb0\xef\xbc\x88list\xef\xbc\x89\npath: \xe5\xb7\xa5\xe4\xbd\x9c\xe8\xb7\xaf\xe5\xbe\x84\xef\xbc\x8c\xe7\x94\xa8\xe4\xba\x8e\xe4\xbf\x9d\xe5\xad\x98\xe4\xb8\xb4\xe6\x97\xb6\xe6\x96\x87\xe4\xbb\xb6\nspecies_name: \xe7\x89\xa9\xe7\xa7\x8d\xe5\x90\x8d\xe7\xa7\xb0\xef\xbc\x8c\xe5\xa6\x82 "Klebsiella pneumoniae"\neg.db: \xe8\xbf\x94\xe5\x9b\x9e\xe5\x8c\x85\xe5\x90\xabgene\xef\xbc\x8cfunction\xe7\x9a\x84pd.DataFrame\n        \n\xe2\x88\xa0( \xe1\x90\x9b \xe3\x80\x8d\xe2\x88\xa0)\xef\xbc\xbf \xe2\x88\xa0( \xe1\x90\x9b \xe3\x80\x8d\xe2\x88\xa0)\xef\xbc\xbf \xe2\x88\xa0( \xe1\x90\x9b \xe3\x80\x8d\xe2\x88\xa0)\xef\xbc\xbf \xe2\x88\xa0( \xe1\x90\x9b \xe3\x80\x8d\xe2\x88\xa0)\xef\xbc\xbf \xe2\x88\xa0( \xe1\x90\x9b \xe3\x80\x8d\xe2\x88\xa0)\xef\xbc\xbf \xe2\x88\xa0( \xe1\x90\x9b \xe3\x80\x8d\xe2\x88\xa0)\xef\xbc\xbf\n\n        '
    print(doc.decode('utf-8'))

# def transposon():
#     tn_url = 'https://transposon.lstmed.ac.uk/tn-registry'
#     with request.urlopen(tn_url) as f:
#         tab = f.read()
#         tab = str(tab)
#     # bb=bs(tab)
#     # bb=bb.find_all('li',"first last")
#     bb=bs(tab)
#     bb=bb.find_all('a',title="Go to last page")
#     last_page = int(str(bb[0]).split('page=')[1].split('"')[0])

#     link = []

#     for i in range(last_page):
#         with request.urlopen(tn_url+"?page="+str(i)) as f:
#             tab = f.read()
#             tab = str(tab)
#         bb=bs(tab)
#         bb=bb.find_all('div', "item-list")
#         link += bb
#         time.sleep(1)

#     print(link)

# transposon()