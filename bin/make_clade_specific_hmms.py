#!/usr/bin/env python
###############################################################################
#
# fundec.py <tree> <database_file> <%id cutoff>
# I should probably include some regular expression recognition to make the
# counts of each annotation more robust. I imagine that as is this
# code will only cluster the major branches. I should really set the
# %ID cutoff high to avoid over clustering. I dont want that to come
# back to bite.
#
# Come to think of it, would some kind of modification of Levenshtein 
# distance work here?
#
###############################################################################
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                             #
###############################################################################
 
__author__ = "Joel Boyd"
__copyright__ = "Copyright 2015"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"
 
###############################################################################

import sys
import logging
import argparse
import os
from run import Hmmm

###############################################################################

debug={1:logging.CRITICAL,
       2:logging.ERROR,
       3:logging.WARNING,
       4:logging.INFO,
       5:logging.DEBUG}

###############################################################################

class CustomHelpFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        return text.splitlines()

    def _get_help_string(self, action):
        h = action.help
        if '%(default)' not in action.help:
            if action.default != '' and \
               action.default != [] and \
               action.default != None \
               and action.default != False:
                if action.default is not argparse.SUPPRESS:
                    defaulting_nargs = [argparse.OPTIONAL,
                                        argparse.ZERO_OR_MORE]

                    if action.option_strings or action.nargs in defaulting_nargs:

                        if '\n' in h:
                            lines = h.splitlines()
                            lines[0] += ' (default: %(default)s)'
                            h = '\n'.join(lines)
                        else:
                            h += ' (default: %(default)s)'
        return h

    def _fill_text(self, text, width, indent):
        return ''.join([indent + line for line in text.splitlines(True)])

def phelp():
    print """
                        ~ Make clade specific HMMs ~

Input options:
    --input_sequences    -    Input sequences to align and use for constructing
                              a tree for partitioning. REQUIRED
    --input_tree         -    Input tree in newick format. NOT REQUIRED

Output options:
    --output_directory   -    Directory to which processing files and results
                              are written to
    --log                -    Output logging information to this file
    --verbosity          -    Level of verbosity, 1 (silent) to 5 (verbose)
    --force              -    Force overwrite a previous run with the same name
    
Runtime options:
    --percentile         -    The percentile cutoff to partition the tree with
                              (The default is the 30th percentile, which will 
                              collapse the major groups in the tree, but more 
                              granular partitioning can be achieved using a 
                              lower percentile)
                              
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                
                                  ABOUT
This tool is used to create clade-specific HMMs from a phylogenetic tree. The
tree can be provided, or created in the program using FastTreeMP. The HMMs are 
taxonomy agnostic, meaning that the only criteria used to group sequences 
together to create a HMM is phylogenetic distance. HMMs therefore represent 
phylogenetic, not taxonomic groupings. HMMs that are created are tested on the 
sequences provided, not an independent dataset, so take the calibration results 
with a grain of salt. That said, they should provide a rough guideline of the 
error rate of each HMM (suggested cutoffs are provided). 

This program does take time to run on large gene families. Be patient.
""" 

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='''Automatically partition a tree''')
    parser.add_argument('--input_tree')
    parser.add_argument('--input_sequences', required=True)
    parser.add_argument('--output_directory',
                        default = "partition_output")
    parser.add_argument('--force',
                        action = 'store_true',
                        default = False)
    parser.add_argument('--percentile', type = int, default = 30)
    parser.add_argument('--log')
    parser.add_argument('--verbosity', type = int, default = 4)

    if(len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help'):
        phelp()
    else:

        args = parser.parse_args()
    
        if args.log:

            logging.basicConfig(filename=args.log, level=debug[args.verbosity], format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
        else:
            logging.basicConfig(level=debug[args.verbosity], format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
        logging.info(' '.join(sys.argv))
        hmmm = Hmmm(args)
        hmmm.run()
        