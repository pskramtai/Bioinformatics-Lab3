import sys
import pandas as pd
import matplotlib.pyplot as plt
import collections
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def find_encoding_scheme():
    with open("reads_for_analysis.fastq", 'r') as f:
        line = f.readline()
        
        min_ascii = None
        max_ascii = None

        while line:
            if line[0] == '+':
            
                line = f.readline()
                line = line.strip()
                
                for c in line:
                    ascii_code = ord(c)
                    
                    if min_ascii is None or ascii_code < min_ascii:
                        min_ascii = ascii_code
                    if max_ascii is None or ascii_code > max_ascii:
                        max_ascii = ascii_code
                
            line = f.readline()
        
        print(min_ascii)
        print(max_ascii)
        if min_ascii >= 33 and max_ascii <= 126 and max_ascii-33 == 40:
            print("Sanger Phred+33")
        elif min_ascii >= 59 and max_ascii <= 126 and max_ascii-64 == 40:
            print("Solexa Solexa+64")
        elif min_ascii >= 64 and max_ascii <= 126 and max_ascii-64 == 40:
            print("Illumina 1.3+ Phred+64")
        elif min_ascii >= 64 and max_ascii <= 126 and max_ascii-64 == 41:
            print("Illumina 1.5+ Phred+64")
        elif min_ascii >= 33 and max_ascii <= 126 and max_ascii-33 == 41:
            print("Illumina 1.8+ Phred+33")
        else:
            print("Encoding scheme not detected")


def get_df():
    with open("reads_for_analysis.fastq", 'r') as f:
        line = f.readline()
        lineId = 0
        
        read = 0
        
        test = []
        while line:
            if line[0] == '+':
                line = f.readline()
                lineId +=1
            elif line[0] == '@':
                line = f.readline()
                lineId +=1
                line = line.strip()

                read += 1
                c_count = line.count('C')
                g_count = line.count('G')

                counter = collections.Counter(line)


                cg_distribution = (c_count + g_count) / len(line)
            


                test.append((lineId,cg_distribution))
                
            line = f.readline()
            lineId+=1

    
        df = pd.DataFrame(test,columns=['id','ratio'])
        return df

def blast_search(df):
    with open("reads_for_analysis.fastq", 'r') as f:
        lines = f.readlines()
        print(df)
        for index, row in df.iterrows():
            lineId = int(row['id'])
            line = lines[lineId]
            line = line.strip()
            print(line)
            result_handle = NCBIWWW.qblast("blastn", "nr", line, entrez_query="bacteria[ORGN]")
            print("test")
            blast_results = NCBIXML.read(result_handle)
            for alignment in blast_results.alignments:
                for hsp in alignment.hsps:
                    print("sequence:", alignment.title)

find_encoding_scheme()
df = get_df()
df.plot.bar(width=1)
plt.show()

df = df.sort_values("ratio", ascending=False)
blast_search(df[0:5])




