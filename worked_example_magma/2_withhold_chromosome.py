import pandas as pd
import sys

trait = sys.argv[1]

filename = '%s.genes.raw' %(trait)

for chrom in range(1,23):
    outname = '%s_noChr%s.genes.raw' %(trait, str(chrom))
    with open( filename) as x:
        with open(outname, 'w') as out:
            first_header = x.readline()
            second_header = x.readline()
            out.write( first_header )
            out.write( second_header )

            for line in x:
                gene_num, listed_chrom, start, stop = line.strip().split( )[0:4]
                if str(listed_chrom) != str(chrom):
                    out.write( line )
