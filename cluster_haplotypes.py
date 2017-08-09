"""
TODO: Find better way of assigning numbers... possibly having to do with their numeric values?
"""


import numpy
import math
import argparse
from collections import Counter
from functools import reduce

class kMeansHaplotypeCluster:
	
	"""
	Initializes k means clustering algorithm
	"""
	def __init__(self):
		self.symbols = None
		self.num_iterations = None
		self.optimize = None
		self.method = None


	"""
	Run algorithm
		raw_haplotypes: iterable with haplotypes (any symbols appropriate)
		max_clusters: maximum number of clusters to attempt
	"""
	def run(self, raw_haplotypes, max_clusters, num_iterations=1000, optimize=True, method=''):
		self.method = method
		self.num_iterations = num_iterations
		self.optimize = optimize

		if self.symbols is None:
			self.symbols = self.__generate_symbol_dict__(raw_haplotypes) # generate symbol dictionary
		haplotypes = [[self.symbols[s] for s in hap] for hap in raw_haplotypes] # translate haplotypes into symbols
		clusters = None
		if optimize:
			clusters = self.__optimize_cluster_distances__(haplotypes, max_clusters)	# run clustering algorithm
		else:
			clusters = self.__k_means__(haplotypes, max_clusters)[0]
		return [[[self.symbols[s] for s in hap] for hap in hap_list] for hap_list in clusters] # retranslate into original symbols and return
	
	"""
	generate symbol dictionary
	"""

	def getSymbolDict(self, raw_haplotypes):
		if self.symbols is None:
			self.symbols = self.__generate_symbol_dict__(raw_haplotypes) # generate symbol dictionary
		return self.symbols
	"""
	Generates consensus sequence from haplotypes
	"""
	def __generate_consensus_sequence__(self, haplotypes):
		consensus = list()
		arr = numpy.array(haplotypes, order='F') #generate array from list of haplotypes
		for col in arr.T:
			consensus.append(Counter(col).most_common(1)[0][0]) #add most common nucleotide in each position to consensus sequence
		return consensus
		

	"""
	Score sequence (sum of squares)
	"""
	def __score_sequence__(self, seq, consensus):
		if not len(seq) == len(consensus):
			raise ValueError('Sequence: %s and consensus: %s lengths do not match' % (seq, consensus))
		return reduce(lambda x, y : x+y, [(seq[i] - consensus[i])**2 for i in range(len(seq))])


	"""
	Generate symbol values for haplotypes
	"""
	def __generate_symbol_dict__(self, haplotypes):
		symbols = dict()
		seen = dict()
		num = 1
		for hap in haplotypes:
			for sym in hap:
				sorted_sym = ''.join(sorted(list(sym)))
				if sorted_sym not in seen:
					symbols[sym] = num
					symbols[num] = sym
					seen[sorted_sym] = [sym]
					num += 1
				elif sym not in symbols:
					symbols[sym] = symbols[seen[sorted_sym][0]]
					seen[sorted_sym].append(sym)
				"""
				if not sym in symbols:
					symbols[sym] = num
					symbols[num] = sym
					num += 1
				"""
		return symbols

	"""
	Average distance from cluster centroid
	"""
	def __avg_cluster_distance__(self, haplotypes, consensus=None):
		if consensus is None:
			consensus = self.__generate_consensus_sequence__(haplotypes)
		return reduce(lambda x,y : x+y, [self.__score_sequence__(seq, consensus) for seq in haplotypes])/len(haplotypes)


	"""
	Use k-means algorithm to cluster haplotypes
	"""
	def __k_means__(self, haplotypes, num_clusters):
		clusters = [[] for i in range(num_clusters)] # list of lists to store clustered records
		#centroids = [haplotypes[i] for i in numpy.random.choice(len(haplotypes), size=num_clusters, replace=False)] #Forgy method of initialization (random choice)
		#use modified Forgy method (unique random choice) for initialization
		
		centroids = []
		tries = 0
		if self.method == 'partition':
			haps_copy = list(haplotypes) #copy list of haplotypes
			while len(haps_copy) > 0:
				index = numpy.random.choice(len(haps_copy), size=1, replace=False)[0] #get random haplotype
				clusters[tries].append(haps_copy.pop(index)) #remove from list and add to appropriate cluster
				# move to next cluster in a circular manner
				tries += 1 
				if tries >= len(clusters):
					tries = 0
			centroids = [self.__generate_consensus_sequence__(cluster) for cluster in clusters]
		elif self.method == 'forgy':
			while tries < self.num_iterations:
				index = numpy.random.choice(len(haplotypes), size=1, replace=False)[0]
				hap = tuple(haplotypes[index])
				if hap in centroids:
					tries += 1
					continue
				centroids.append(hap)
				if len(centroids) == num_clusters:
					break
		else:
			centroids = self.__k_means_plus_plus__(haplotypes, num_clusters)

		if len(centroids) < num_clusters:
			raise ValueError('Could not determine %d distinct centroids for starting conditions' % num_clusters)
		
		for i in range(self.num_iterations):
			clusters = [[] for i in range(num_clusters)] # list of lists to store clustered records
			for hap in haplotypes: # iterate over haplotypes
				distances = enumerate([self.__score_sequence__(hap, centroid) for centroid in centroids]) #calculate distances from each of the centroids
				cluster = min(distances, key=lambda x:x[1]) # get the minimum distance
				clusters[cluster[0]].append(hap) # append to the cluster with the minimum distance
			for index, cluster in enumerate(clusters):
				if cluster:
					centroids[index] = self.__generate_consensus_sequence__(cluster)
		return clusters, centroids


	def __k_means_plus_plus__(self, haplotypes, num_clusters):
		#randomly choose the first centroid
		centroids = []
		index = numpy.random.choice(len(haplotypes), size=1, replace=False)[0]
		centroids.append(haplotypes[index])
		#get the rest of the centroids
		for i in range(num_clusters-1):
			probs = list()
			for hap in haplotypes:
				distances = [self.__score_sequence__(hap, centroid) for centroid in centroids] #calculate distances from each of the centroids
				probs.append(min(distances)) # get the minimum distance and append to the probabilities
			total = sum(probs)
			if total == 0:
				break
			probs = list(map(lambda x : x/total, probs)) #generate probability distribution
			#print(probs)
			
			index = numpy.random.choice(len(haplotypes), size=1, replace=False, p=probs)[0]
			centroids.append(haplotypes[index])
		#print(centroids)
		return centroids
	

	"""
	Find the minimum average cluster distance for 2 clusters to max_clusters
	"""
	def __optimize_cluster_distances__(self, haplotypes, max_clusters):
		min_avg_distance_clusters = None
		min_avg_distance = math.inf
		for max_c in range(2, max_clusters + 1):
			try:
				clusters, centroids = self.__k_means__(haplotypes, max_c)
			except ValueError:
				break
			distance = reduce(lambda x,y : x+y, [self.__avg_cluster_distance__(cluster, consensus=centroid) for cluster, centroid in zip(clusters, centroids) if cluster])/max_c #get average of average distances between clusters
			if distance < min_avg_distance:
				min_avg_distance_clusters = clusters
				min_avg_distance = distance
		return min_avg_distance_clusters #return cluster with minimum average distance from cluster centroid


def __main__():
	SEP = ' ' #separator
	parser = argparse.ArgumentParser('Cluster Haplotypes')
	parser.add_argument('-o', '--output', help='Output file for clusters')
	parser.add_argument('-c', '--clusters', help='Maximum number of clusters', type=int, required=True)
	parser.add_argument('-n', '--iterations', help='Number of iterations', type=int, default=1000)
	parser.add_argument('-m', '--method', help='Method to use (default kmeans++)', type=str)
	parser.add_argument('--nooptimize', help='Do not optimize cluster number', action='store_true')
	parser.add_argument('input', help='Input filename')
	args = parser.parse_args()
	optimize = not args.nooptimize
	clust = kMeansHaplotypeCluster()
	
	#haplotypes = ['ATCGGTGCACAATGG','ATGGGTGCACAAGGG','ATTGGTGCACAATGG','ATGGGTGCACAATGG','ATGGGTGCACAAGGG','ATCGGTGCACAATGG','ATGGGTGCACAATGG','ATCGGTGCACAATGG','ATGGGTGCACAAGGG','ATCGGTGCACAATGG','ATCGGTGCACAAGGG']
	#print(clust.run(haplotypes, 3))
	
	haplotypes = list()
	header = ''
	samples = dict()
	temp_sample_list = list()
	with open(args.input, 'r') as f:
		for line in f:
			line = line.strip()
			if line[0] == '#':
				header = line
			else:
				data = line.split('\t')
				hap = data[0].split(SEP)
				#hap = list(data[0])
				haplotypes.append(hap)
				temp_sample_list.append(data[1])
	symbol_dict = clust.getSymbolDict(haplotypes)
	haplotypes = [[symbol_dict[symbol_dict[sym]] for sym in hap] for hap in haplotypes]

	for hap, hap_samples in zip(haplotypes, temp_sample_list):
		samples[' '.join(hap)] = hap_samples
	
	del temp_sample_list
	
	if args.output:
		with open(args.output, 'w') as wf:
			wf.write(header + '\n')
			for cluster in clust.run(haplotypes, args.clusters, num_iterations=args.iterations, optimize=optimize, method=args.method): 
				wf.write('\n'.join([' '.join(hap) + '\t' + samples[' '.join(hap)] for hap in cluster]) + '\n')
				wf.write('\n\n')
	else:
		print(header)
		for cluster in clust.run(haplotypes, args.clusters, num_iterations=args.iterations, optimize=optimize, method=args.method):
			print('\n'.join([' '.join(hap) + '\t' + samples[' '.join(hap)] for hap in cluster]))
			print('\n')


if __name__ == '__main__':
	__main__()
