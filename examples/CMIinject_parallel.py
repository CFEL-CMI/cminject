#!/opt/local/bin/python
#simple multiprocessing wrapper for CMIinject
# Daniel Horke, 2018

# import generally necessary modules
import subprocess
import numpy as np
from multiprocessing import Pool
import math
import os
import sys
import argparse
from time import strftime

parser = argparse.ArgumentParser()
parser.add_argument(
                    "-f",
                    "--flowfield",
                    type=str,
                    help="Flowfield to use for simulations")
parser.add_argument(
                    "-d",
                    "--directory",
                    type=str,
                    default = 'test',
                    help="Output directory. default = test ")
parser.add_argument(
                    "-v",
                    "--velocity",
                    type=float,
                    default=-100.0,
                    help="Specify center of z velocity distribution. Default is -100.00")
parser.add_argument(
                    "-p",
                    "--particles",
                    type=int,
                    default=1000,
                    help="Specify number of particles per job. Default is 1000.")
parser.add_argument(
                    "-c",
                    "--cores",
                    type=int,
                    default=2,
                    help="specify number of parallel processes (=no of available cores). default = 2")
parser.add_argument(
                    "-j",
                    "--jobs",
                    type=int,
                    default=2,
                    help="Number of jobs to run. Total number of particles flow = particles*jobs. default = 2")
args=parser.parse_args()

sys.stderr = open(os.devnull, 'w')

print(strftime("%Y-%m-%d %H:%M:%S"), ": Simulating a total of "+str(args.jobs*args.particles)+"  particles on "+str(args.cores)+" nodes.")

#define function call for flying particles. This will be called args.jobs number of times
def worker(i):
    print(strftime("%Y-%m-%d %H:%M:%S"), ": Running job "+str(i+1)+"  of "+str(args.jobs)+"!")
    subprocess.call("python run.py "+str(args.flowfield) +" "+str(args.velocity)+" "+str(args.directory)+"/"+str(i).zfill(4)+"_ "+str(args.particles), shell=True)

def main(argp=None):
    # Create the number of threads you want
    pool = Pool(args.cores);
    try:
        # Spawn up to 9999999 jobs, I think this is the maximum possible.
        pool.map_async(worker, (i for i in range(args.jobs))).get(9999999)
    except KeyboardInterrupt:
        print('User interrupt.')
    pool.close()

if __name__ == '__main__':
    main()
