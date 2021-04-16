# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 14:14:40 2018

@author: s1767711

"""


map_data=open("Cervus_elaphus_SNP-linkage.ordered_qced.txt")
new_SNP_names=open("Cervus_elaphus_SNP_names_with_linkage.txt", "w")
for line in map_data:
     line=line.rstrip("\n")
     line_list=line.split("\t")
     snp=line_list[0]
     link=line_list[1]
     new_SNP_names.write(snp+"-"+link+"\n")
new_SNP_names.close()


"""
create files containing SNPs of 20 SNP windows along each autosome
"""

import re
SNPs=open("Cervus_elaphus_SNP_names_with_linkage.txt").readlines()
for n in range(1,34):
    list_link=list()
    pattern="-"+str(n)+"$"
#    print(pattern)
    for line in SNPs:
        line=line.rstrip("\n") 
        match=re.search(pattern, line)
        if match:
           list_link.append(line) 
#    print(len(list_link))
    window_size=20
    step=10
    range_stop=step-1
    positions=list(range(0, len(list_link)-range_stop, step))
    for pos in positions:
        window=list_link[pos:pos+window_size]
        if len(window)<20:
            window=list_link[len(list_link)-21:-1]
        f= open("window_link"+str(n)+"_pos"+str(pos)+".txt",'w') 
        for snp in window:
            snp_split=snp.split("-")
            snp=snp_split[0]
            f.write(str(snp)+"\n")
        f.close()
        
        
 