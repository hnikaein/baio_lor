import os
import subprocess
from multiprocessing.pool import Pool

FNULL = open(os.devnull, 'w')
ref_name = "grch38.fasta"
refs_count = 197163
log_chunk = 15
processes = 32


def run_aryana(i):
    subprocess.call(["../aryana/aryana", "index", "chunks/%s_%d_%d.fasta" % (ref_name, log_chunk, i)], stderr=FNULL)
    subprocess.call(["../aryana/aryana", "fa2bin", "chunks/%s_%d_%d.fasta" % (ref_name, log_chunk, i)], stderr=FNULL)
    if i % 1000 == 0:
        print("------->indexing %06d" % i)


p = Pool(processes)
p.map(run_aryana, range(refs_count))
p.close()
