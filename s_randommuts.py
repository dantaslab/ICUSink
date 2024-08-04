#!/bin/env python3

## Author: JooHee Choi, joohee.choi@wustl.edu, used by Erin Newcomer for ICU Sink project.
##this script will generate a random number of mutations throughout the genome and note for each gene, if is mutated X times 
## Usage: python3 s_randommuts.py
##output file: randommuts.csv 

import random


gff=[]
strains=[]
sizes={} #dictionary of strain: genomesize
allstarts={} #dictionary of strain: [starts]
mutdict={} #dictionary of strain: iso count 
matches={} #library of patient: reference gff file 

def main():

	snp_m=open('patients.txt', 'r') # list of patients (matrix) names
	gff_f=open('refgff.txt', 'r') # list of reference gff files _.gff
	
	for line in gff_f:
		gffname=line.strip()
		gff.append(gffname) 

	gff_f.close()

	for line in snp_m:
		strain=line.strip()
		strains.append(strain)

	snp_m.close()
	
	n=len(gff)

	##set up corresponding matches 
	for i in range(n):
		matches[strains[i]]=gff[i]
	


	for strain in strains: 
		starts, genomesize=strainlength(strain)
		sizes[strain]=genomesize
		allstarts[strain]=starts
		mutdict[strain]=mutcount(strain)

	#print(mutdict)




def strainlength(strain):  ##calculate length of genomes for each strain 

	ref=matches[strain]
	file=open('gff/'+ref)
	nodestarts=[] #list of starting positions for each node 
	counter=0
	nodestarts.append(1) ## to start the list 
	#nodecounter=0

	for line in file:
		if line[0:3]=='##s':
			line_list=line.strip().split(' ')
			counter+=int(line_list[3])
			nodestarts.append(counter)
			#nodecounter+=1


	nodestarts=nodestarts[:-1] ##get rid of last entry, which is genome size

	file.close()
	return(nodestarts, counter)


###function to get number of mutations per strain 
def mutcount(strain):
	
	
	
	matrix_file='snp_matrix/'+ strain+'_merged_parsed.txt' ##snp count ###CHANGE THIS LINE


	mfile=open(matrix_file, 'r')
	samplist=[]
	

	first='T'
	for line in mfile:
		if first=='T':
			if line[0]!='/': ## excluding for sample names 

				sampleline=line.strip().split('\t')
				n=len(sampleline)

				first='F'
			else:
				n=0

		if first=='F':

			break

	mfile.close()

	

	##added on second portion to add indels 
	matrix_file2='indel_matrix/'+strain+'_indel_merged_parsed.txt' ##indel count 
	
	mfile=open(matrix_file2, 'r')
	samplist=[]
	

	first='T'
	for line in mfile:
		if first=='T':
			if line[0]!='/':

				sampleline=line.strip().split('\t')
				m=len(sampleline)
				first='F'

			else:
				m=0
				

		if first=='F':
			break



	mfile.close()
	
	count=n+m
	
	return count

	


def matchgff(strain): ##create library of genes at each node 

	ref=matches[strain]
	file=open('gff/'+ref)
	gff_lib={} ##node:[[start,end,gene],[start,end,gene],....,] 
	
	for line in file:
		
		if (line[0]!='#' and line[0]!='>' and line[0]!='A' and line[0]!='T' and line[0]!='G' and line[0]!='C'):
			line_list=line.strip().split('\t')

				
			info=str(line_list[8])

			if 'gene=' in info:
				temp=info.split('gene=')[1]
				gene=temp.split(';')[0]

			else:
				gene='hypothetical protein'

			temp=[]
			node=line_list[0] ##so that node name is NODE_70 

				
			temp.append(line_list[3]) #start
			temp.append(line_list[4]) #end 
			temp.append(gene) ##gene.location of 'gene' in info tab too inconsistent, more over some don't even have it
		
			if node not in gff_lib:
				gff_lib[node]=[]
				gff_lib[node].append(temp)
			else:
				gff_lib[node].append(temp)

	file.close()

	return gff_lib



##match between node, nodepos and resulting gene 
def matchnode(nodenum, pos):
	nodename=str(nodenum)

	if nodename in gfflib: 
			
		list_pos=gfflib[nodename]  ##list of [start, end, gene], [start, end, gene], 
		#print(list_pos)

		m=len(list_pos)
		gene=''
		n=0
		for k in range(m):
			
			if int(list_pos[k][0]) < pos and int(list_pos[k][1]) > pos: # any hit 
				n+=1

				if n==1: ##first hit 
					gene+=str(list_pos[k][2])

				else:
					gene=gene+','+str(list_pos[k][2])

		'''
		if gene=='': #because no matches were found at all

			temp='intron'
		#		gene='intron'	
			for j in range(1000):
				pos2=pos+j
				pos3=pos-j
					

				if int(list_pos[k][0]) < pos2 and int(list_pos[k][1]) > pos2: # any hit 
					n+=1
						#m+=1
					if n==1: ##first hit 
						temp=temp+'_ds'+str(j)+'_'+str(list_pos[k][2])
						break

								
				if int(list_pos[k][0]) < pos3 and int(list_pos[k][1]) > pos3: # any hit 
					n+=1

					if n==1: ##first hit 
						temp=temp+'_us'+str(j)+'_'+str(list_pos[k][2])
						break


			gene=temp
		'''
		if gene=='':
			gene='intron'
	else:
		gene='unknown node'
	

	return gene




if __name__ == '__main__':
    main()




n=len(strains) #patients 




shared_sum=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]##list for keeping track of how many genes are shared by how many strains 
s=len(shared_sum)

itera=1000
count_iter=itera

all_muts_sum={}
outfile=open('randommuts.csv','w') ##write a new file 

while count_iter > 0:
	all_muts={} ##dictionary of gene: how many strains have that mutation (not counting intronic for now)


	for i in range(n):
		
		strain_muts=[] ##list of genes mutated in this strain , not always == nummuts since there may be nonunique entries or SNPs with multiple genes 

		strain=strains[i]

		genomesize=sizes[strain]
		starts=allstarts[strain]


		gfflib=matchgff(strain) #get gff lib 

		nodecount=len(starts)



		mutnum=mutdict[strain] #total number of mutations to generate. placement here important since while loop 


		##generate a certain number of random numbers 
		while mutnum >0:
			query=random.randrange(genomesize)
			mutnum-=1

			nodenum=0
			nodepos=0

			##take our number and figure out the node and position on that node
			for k in range(nodecount-1):
				if query>= starts[k] and query < starts[k+1]:
					nodenum=k+1
					nodepos=query-starts[k]+1



			if query >= starts[nodecount-1]: ##if belongs to LAST segment of node 
				nodenum=nodecount
				nodepos=query-starts[nodecount-1]+1

			#print(query, nodenum, nodepos, genomesize)
			##map back this node num and pos into gene info 
			match=matchnode(nodenum, nodepos)

			##in case match is a list of genes appended by ','

			genelist=match.split(',')
			for gene in genelist:
					if gene not in strain_muts: ##excluding introns for now 
						strain_muts.append(gene)

		
		##add the genes into strain_muts dictionary 
		for mut in strain_muts:
			if mut in all_muts:
				all_muts[mut]+=1
			else:
				all_muts[mut]=1
		#print(strain, strain_muts)


	shared=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] ##list for keeping track of how many genes are shared by how many strains 
	



	for entry in all_muts:
		if entry!='intron' and entry!='hypothetical protein' and entry!='unknown node':
			k=all_muts[entry] ##k is how many strains share that gene mutation
			shared[k-1]+=1 ## adds 1 to the position in 'shared' corresponding to that number of strains (ex: shared by 1 strain, pos=0)

			outfile.write(str(count_iter))
			outfile.write(',')
			outfile.write(entry)
			outfile.write(',')
			outfile.write(str(k))
			outfile.write('\n')
	

	for i in range(s):
		shared_sum[i]+=shared[i]


	print(count_iter)
#	print(all_muts)
	print(shared)


	for j in all_muts: ##counting mutated strains per gene 
		if j in all_muts_sum:
			all_muts_sum[j]+=all_muts[j]
		else:
			all_muts_sum[j]=all_muts[j]

	#print(all_muts_sum)
	count_iter-=1

outfile.close()

for i in range(s):
	shared_sum[i]=shared_sum[i]/itera

print(shared_sum)

for j in all_muts_sum:
	print(j, all_muts_sum[j]/itera, sep=',')
