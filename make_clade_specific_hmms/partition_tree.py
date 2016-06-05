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
# System imports
import logging
import numpy as np
import random

from itertools import combinations
from skbio.tree import TreeNode

###############################################################################
############################### - Exceptions - ################################

class MalformedTreeException(Exception):
    pass

###############################################################################
################################ - Classes - ##################################

class Cluster:
    
    UNCLUSTERED = "UNCLUSTERED"
    
    def _node_dist(self, node):
        '''
        Returns a list of distances tip-to-tip from the node provided. 
        
        Parameters
        ----------
        node: skbio TreeNode obj
            http://scikit-bio.org/docs/latest/generated/skbio.tree.TreeNode.html#skbio.tree.TreeNode
        Returns
        -------
        Array of distances from tip-to-tip
        '''
        logging.debug("Calculating tip-to-tip distances for node: %s" \
                                                                % node.name)
        distances=[]
        if len(list(node.tips())) > 200:
            
            for tip_a, tip_b in list(random.sample(list(combinations(node.tips(), 2)), 6000)):
                distances.append(tip_a.distance(tip_b))
            distances=np.array(sorted(distances))
        else:
            for tip_a, tip_b in combinations(node.tips(), 2):
                distances.append(tip_a.distance(tip_b))
            distances=np.array(sorted(distances))
        return distances
    
    def _rename(self, node, name):
        
        if node.name:
            try: 
                float(node.name)
                node.name = "%s:%s" % (node.name,
                                       name)
            except:
                node.name = "%s; %s" % (node.name,
                                        name)
        else:
            node.name = name
                   
    def depth_partition(self, input_tree, percentile, output_tree):
        '''
        Attempt to cluster tree with nodes of tip-to-tip distrubution <
        an nth percentile cutoff of the whole-tree distance distribution. 
        A better description can be found in the citation below.
        
        Parameters
        ----------
        tree: skbio TreeNode obj
            http://scikit-bio.org/docs/latest/generated/skbio.tree.TreeNode.html #skbio.tree.TreeNode
                
        percentile: float
            The percentile cutoff to use to determine the cutoff from clading
            from a given node.
        
        Clustering method modified from Prosperi et al method:
        Prosperi, M.C.F., et al. A novel methodology for large-scale phylogeny 
        partition. Nat. Commun. 2:321 doi: 10.1038/ncomms1325 (2011).
        
        http://www.nature.com/ncomms/journal/v2/n5/full/ncomms1325.html
        '''
        tree = TreeNode.read(input_tree)
        
        cluster_count = 1
        clustered = set()
        clusters = {}
        logging.debug("Calculating %ith percentile cutoff from root" \
                                                % (percentile))
        whole_tree_distribution = self._node_dist(tree)
        
        cutoff = np.percentile(whole_tree_distribution, percentile)
        logging.debug("Cutoff (%ith percentile): %f" % (percentile,
                                                        cutoff))
        for node in tree.preorder():
            if node in clustered:
                continue
            elif node.is_tip():
                continue
            else:
                node_distribution = self._node_dist(node)
                median=np.median(node_distribution)
                logging.debug("Median of node: %f" % median)
                if median <= cutoff:
                    logging.debug("Cluster found!")
                    cluster_name =  "partition_%i" % (cluster_count)
                    clusters[cluster_name] = [x.name.replace(' ','_') 
                                              for x in node.tips()]
                    self._rename(node, cluster_name)
                    cluster_count+=1
                    for descenent in node.traverse():
                        clustered.add(descenent)
        logging.info("%i depth cluster(s) found in tree" % (cluster_count-1))
        tree.write(output_tree, "newick")
        
        logging.debug("Recording tips that were not partitioned")
        clusters[self.UNCLUSTERED] = []
        for tip in tree.tips():
            if tip not in clustered:
                clusters[self.UNCLUSTERED].append(tip.name.replace(' ','_'))
        return clusters
###############################################################################
###############################################################################
###############################################################################
###############################################################################
