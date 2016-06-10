#!/usr/bin/env python
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

import os
import tempfile
import subprocess
import shutil
import logging
from partition_tree import Cluster
from collections import Counter
from Bio import SeqIO

###############################################################################
###############################################################################
###############################################################################
###############################################################################
    
class Hmmm: 
    
    ALIGNMENTS_SUBDIRECTORY = "alignments"
    HMMS_SUBDIRECTORY = "hmms"
    
    def __init__(self, args):
        logging.info("Preparing run")
        self.percentile = args.percentile

        if os.path.isdir(args.output_directory):
            if args.force:
                logging.warning("Deleting directory: %s" % args.output_directory)
                shutil.rmtree(args.output_directory)
            else:
                raise Exception("Directory %s already exists!" \
                                        % args.output_directory)

        self.output_directory = args.output_directory
        self.alignments_subdirectory = os.path.join(self.output_directory, 
                                                    self.ALIGNMENTS_SUBDIRECTORY)
        self.hmms_subdirectory = os.path.join(self.output_directory, 
                                              self.HMMS_SUBDIRECTORY)
        os.mkdir(self.output_directory)
        os.mkdir(self.alignments_subdirectory)
        os.mkdir(self.hmms_subdirectory)
        
        self.sequences = args.input_sequences
        self.base = os.path.splitext(os.path.basename(self.sequences))[0]
        
        if not args.input_tree:        
            logging.info("No tree was provided. Alinging sequences (mafft), then constructing a phylogenetic tree (FastTreeMP) on the fly")      
            output_alignment = os.path.join(self.output_directory,
                                            "%s.aln.fasta" % self.base)
            output_tree      =  os.path.join(self.output_directory,
                                             "%s.nwk" % self.base)
            self._make_align(self.sequences, output_alignment)
            self._make_tree(output_alignment, output_tree)
            self.tree = output_tree
        else:              
            self.tree = args.input_tree   
        
        shutil.copy(self.sequences, 
                    os.path.join(self.output_directory,
                                 self.sequences))
        self.sequences = os.path.join(self.output_directory, self.sequences)
        self.sequences_dict = SeqIO.to_dict(SeqIO.parse(self.sequences, "fasta"))
        self.output_tree = os.path.join(self.output_directory,
                                        "%s.partitioned%s" \
                                            % os.path.splitext(os.path.basename(self.tree)))
    
    def _make_tree(self, sequences, output_tree):
        logging.debug("Contstructing phylogenetic tree using FastTreeMP")
        cmd = "FastTreeMP -quiet %s > %s 2>/dev/null" % (sequences, output_tree)
        logging.debug("Running: %s" % cmd)
        subprocess.call(cmd, shell=True)
    
    def _make_align(self, name, output):
        
        cmd = "mafft --thread 5 --quiet --anysymbol %s > %s" \
                                % (name, output)
        logging.debug("Running: %s" % cmd)
        subprocess.check_call(cmd, shell=True)
    
    def _align_sequences(self, sequence_group, sequence_names):

        logging.debug("Aligning %i sequences in cluster: %s" \
                                    % (len(sequence_names), sequence_group))
        output_alignment_path = os.path.join(self.alignments_subdirectory,
                                             "%s.aln.fasta" % sequence_group)
        sequence_seqs = [self.sequences_dict[name] for name in sequence_names]
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            SeqIO.write(sequence_seqs, fasta, "fasta")
            fasta.flush()
            self._make_align(fasta.name, output_alignment_path)
        return output_alignment_path
        
    def _create_hmm(self, sequence_group, alignment_file):
        logging.debug("Generating clade-specific HMM for %s" % sequence_group)
        output_hmm_path = os.path.join(self.hmms_subdirectory,
                                       "%s.hmm" % sequence_group)
        cmd = "hmmbuild -o /dev/null %s %s" % (output_hmm_path, alignment_file)
        subprocess.check_call(cmd, shell=True)
        
        return output_hmm_path

    def _check_hmmsearch_output(self, hmmsearch_output):
        result_dictionary = {}
        hmmsearch_output_list = hmmsearch_output.strip()\
                                                .split('\n')
        
        for line in hmmsearch_output_list:
            if line.startswith('#'): continue
            split_line = line.split()
            
            sequence_name = split_line[0]
            hmm_name = os.path.splitext(split_line[3])[0]
            bit_score = float(split_line[7])
            
            if sequence_name in result_dictionary:
                if bit_score > result_dictionary[sequence_name][0]:
                    result_dictionary[sequence_name] = [bit_score, 
                                                        hmm_name]
            else:
                result_dictionary[sequence_name] = [bit_score, 
                                                    hmm_name]
        
        return result_dictionary
    
    def _interpret_results(self, 
                           whole_sequence_results, 
                           fragmented_sequence_results, 
                           number_of_fragments_dict,
                           partitions,
                           whole_output_file,
                           fragmented_output_file):
        
        EXPECTED = "expected"
        OBSERVED = "observed"
        FALSE = "false"
        NOISE_CUTOFF = "noise_cutoff"
        TRUSTED_CUTOFF = "trusted_cutoff"
        
        whole_interpreted_results = {partition_name: {EXPECTED:len(partition_list), 
                                                      OBSERVED:0,
                                                      FALSE: [],
                                                      NOISE_CUTOFF:0,
                                                      TRUSTED_CUTOFF:0} 
                               for partition_name, partition_list in partitions.items()}

        fragmented_interpreted_results = {partition_name: {EXPECTED:len(partition_list), 
                                                      OBSERVED:0,
                                                      FALSE: [],
                                                      NOISE_CUTOFF:0,
                                                      TRUSTED_CUTOFF:0} 
                               for partition_name, partition_list in partitions.items()}

        reversed_partitions = {}        
        for partition_name, sequence_names_list in partitions.items():
            for name in sequence_names_list:
                    reversed_partitions[name] = partition_name
        
        for sequence_name, result in whole_sequence_results.items():
            observed_partition = '_'.join(result[1].split('_')[-2:])
            expected_partition = reversed_partitions[sequence_name]
            if expected_partition == Cluster.UNCLUSTERED: continue
            if expected_partition == observed_partition:
                whole_interpreted_results[observed_partition][OBSERVED]+=1
                if result[0] > whole_interpreted_results[observed_partition][TRUSTED_CUTOFF]:
                    whole_interpreted_results[observed_partition][TRUSTED_CUTOFF]=result[0] 
            else:
                whole_interpreted_results[expected_partition][FALSE].append(
                                                                            '_'.join(result[1].split('_')[1:])
                                                                            )
                if result[0] > whole_interpreted_results[expected_partition][NOISE_CUTOFF]:
                    whole_interpreted_results[expected_partition][NOISE_CUTOFF] = result[0]

        for sequence_name, result in fragmented_sequence_results.items():
            observed_partition = '_'.join(result[1].split('_')[-2:])
            expected_partition = reversed_partitions['_'.join(sequence_name.split('_')[:-1])]
            if expected_partition == Cluster.UNCLUSTERED: continue
            if expected_partition == observed_partition:
                fragmented_interpreted_results[observed_partition][OBSERVED]+=1
                if result[0] > fragmented_interpreted_results[observed_partition][TRUSTED_CUTOFF]:
                    fragmented_interpreted_results[observed_partition][TRUSTED_CUTOFF]=result[0] 
            else:
                fragmented_interpreted_results[expected_partition][FALSE].append(
                                                                                 '_'.join(result[1].split('_')[1:])
                                                                                 )
                if result[0] > fragmented_interpreted_results[expected_partition][NOISE_CUTOFF]:
                    fragmented_interpreted_results[expected_partition][NOISE_CUTOFF] = result[0]
                    
        with open(whole_output_file, 'w') as out_io:
            headers = ["partition_name",
                       "NC",
                       "TC",
                       "recall",
                       "crossover"]
            out_io.write('\t'.join(headers) + '\n')
            
            for partition_name in partitions.keys():
                if partition_name == Cluster.UNCLUSTERED:continue
                results = whole_interpreted_results[partition_name]
                output_line_list = [partition_name,
                                    str(results[NOISE_CUTOFF]),
                                    str(results[TRUSTED_CUTOFF]),
                                    str(round(float(results[OBSERVED])/\
                                        float(results[EXPECTED]),
                                              1))]
                false_result = ''
                for partition, count in Counter(results[FALSE]).items():
                    fraction = str(float(count)/ float(results[EXPECTED]))
                    false_result += "%s,%s" % (partition, fraction)
                output_line_list.append(false_result)   
                out_io.write('\t'.join(output_line_list) + '\n')
        
        with open(fragmented_output_file, 'w') as out_io:
            headers = ["partition_name",
                       "NC",
                       "TC",
                       "recall",
                       "crossover"]
            out_io.write('\t'.join(headers) + '\n')
            
            for partition_name, partition_list in partitions.items():
                if partition_name == Cluster.UNCLUSTERED:continue
                results = fragmented_interpreted_results[partition_name]
                
                total_seqs = 0
                for name in partition_list:
                    total_seqs += number_of_fragments_dict[name]
                
                output_line_list = [partition_name,
                                    str(results[NOISE_CUTOFF]),
                                    str(results[TRUSTED_CUTOFF]),
                                    str(round(float(results[OBSERVED])/\
                                        float(total_seqs),
                                              1))]
                false_result = ''
                for partition, count in Counter(results[FALSE]).items():
                    fraction = str(float(count)/ float(total_seqs))
                    false_result += "%s,%s" % (partition, fraction)
                output_line_list.append(false_result)   
                out_io.write('\t'.join(output_line_list) + '\n')
                
    def _test_hmms(self, hmm_list):
        # TODO: move this to a function of its own
        
        with tempfile.NamedTemporaryFile(suffix='.hmm') as search_hmm:
            for hmm in hmm_list:
                for line in open(hmm):
                    search_hmm.write(line)
                search_hmm.flush()
            logging.debug("Testing recall of whole sequences")
            cmd = "hmmsearch -o /dev/null --domtblout /dev/stdout %s %s" \
                                            % (search_hmm.name, 
                                               self.sequences)
            hmmsearch_output = subprocess.check_output(cmd, shell = True)
            whole_sequence_result = self._check_hmmsearch_output(hmmsearch_output)
            logging.debug("Testing recall of 33 amino acid sequences (equivalent of 100 bp nucletotide read)")
            number_of_fragments_dict = {}
            with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
                for record in self.sequences_dict.values():
                    start_index = 0
                    finish_index = 33
                    sequence_count = 1
                    while finish_index < len(record.seq):
                        fasta.write(">%s_%i\n%s\n" % (record.name,
                                                      sequence_count,
                                                      record.seq[start_index:finish_index]))
                        start_index+=10
                        finish_index+=10
                        sequence_count+=1
                    
                    number_of_fragments_dict[record.name] = sequence_count
                            
                    fasta.flush()
                cmd = "hmmsearch -o /dev/null --domtblout /dev/stdout %s %s" \
                                                % (search_hmm.name, 
                                                   fasta.name)
                hmmsearch_output = subprocess.check_output(cmd, 
                                                           shell = True)
            fragmented_sequence_result = self._check_hmmsearch_output(hmmsearch_output) 
        
        return whole_sequence_result, fragmented_sequence_result, number_of_fragments_dict
        
        
    def run(self):
        c = Cluster()
        hmm_list = []
        
        logging.info("Partitioning tree")
        partitions = c.depth_partition(self.tree,
                                       self.percentile,
                                       self.output_tree)
        logging.info("Aligning clusters and generating clade-specific HMMs")
        for partition_name, partition_sequences in partitions.items():
            if partition_name == c.UNCLUSTERED: continue
            partition_name = "%s_%s" % (self.base, partition_name)
            partition_alignment = \
                    self._align_sequences(partition_name, partition_sequences)
            hmm_file = \
                    self._create_hmm(partition_name, partition_alignment)
            hmm_list.append(hmm_file)
        logging.info("Testing recall of clade-specific HMMs")
        
        whole_sequence_result, fragmented_sequence_result, \
        number_of_fragments_dict = self._test_hmms(hmm_list)
        
        self._interpret_results(whole_sequence_result, 
                                fragmented_sequence_result, 
                                number_of_fragments_dict,
                                partitions, 
                                "%s_whole_sequence_test_results.txt" % os.path.join(self.output_directory,
                                                                     self.base),
                                "%s_fragmented_sequence_test_results.txt" % os.path.join(self.output_directory,
                                                                     self.base))   
        logging.info("Finished")
        