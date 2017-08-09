This is an explanation of the Pipeline used to interpret the G6PD data.
Caleb Palagyi (2017)

1. Input the female data file and male data file into Male2FemaleRankedNC.py
	- This intakes a tab delimited file.
	- Run this program until the graph flatlines or alternates between two points
	- This does not create new genotypes
	- Have to run No Creation otherwise it will take the most common one and pair it with a low result, then it works its way to a singularity.
2. Input the Male results from step 1 and the original female data file into Male2FemaleRanked.py
	- Make sure there is no 000000… result as this will cause every sample to create a new genotype.
	- We removed this sample and then added it back in to avoid this problem.
	- Ranking formula: The sum of the two scores times the Standard Deviation of a Beta Distribution.
3. Took the samples output and clustered them using John Gibson’s algorithm.
	- Outputs the clusters
4. Took the clustered samples and cut the second field using Terminal.
5. Took the file with the second field and used VIM to give each sample it’s own line.
	- Remove all spaces!
6. Created a new file for each cluster using copy, vim, paste.
7. Using the Grepper.py, find the population and gender data for each sample in every cluster.
	- May have to go through the code to select which columns you want
	- This spit out a file consisting of Sample ID, Population, Super Population, Gender
	- Note: You cannot use grep because it will not spit out the data for duplicates.
8. Using the Bash Script Searcher.sh, I counted the occurrences of each population code.
9. Input the results from Searcher.sh into Numbers (Mac Excel) to generate the pie charts.