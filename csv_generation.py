#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse


#defines arguments in the command line
parser = argparse.ArgumentParser()
parser.add_argument("--bed", action="store");
args = parser.parse_args();

# Gets the arguments in the command line
bed_file = args.bed

df = pd.read_csv(bed_file,
                sep = "\t",
                header = None)




