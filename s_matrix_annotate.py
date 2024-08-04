#!/bin/env python3

##Author: JooHee Choi, joohee.choi@wustl.edu, used by Erin Newcomer for ICU Sink project.
##this script will take SNP matrices and write new matrices with gene names. It will also output each mutation per strain.
##usage: python3 s_matrix_annotate.py > output_file_name



strain_dict={} ##'fixed' matrix names: gff file 
strain_matrix=[] #'fixed' matrix names

matrices=[] #list containing all matrix file names  ex) DK-001_1
gff=[] #list containing gff file names ex) xyz.gff

matrices=[] #list of  patient matrices   ex) DK-001_1



def main():
	snp_m=open('patients.txt', 'r') # list of _merged_parsed_txt file names 
	gff_f=open('refgff.txt', 'r') # list of reference gff files _.gff

	
	for line in snp_m:
		matrices.append(line.strip()) ##??   the numbers will get rid of '_merged_parsed_txt'

	
	for line in gff_f:
		gff.append(line.strip()) 

	snp_m.close()
	gff_f.close()


	n =len(matrices)

	for j in range(n):
	
		temp=matrices[j] ###patient name 

		strain_dict[temp]=gff[j]
		strain_matrix.append(temp)




###function to match node list to gff file 
def matchnodes(patient):
	
	########get list of nodes 

	
	matrix_file='snp_matrix/'+ patient+'_merged_parsed.txt'


	mfile=open(matrix_file, 'r')
	nodelist=[]
	parsednodes=[]

	first='T'
	for line in mfile:
		
		if first=='T':
			nodelist=line.strip().split('\t')
			for i in nodelist:
					node=i.split('_')[0]
					pos=i.split('_')[1]
					code=node+'@'+pos ##add in something that's easy to parse
					parsednodes.append(code)

		first='F'


	mfile.close()
	
	return parsednodes


def matchgff(name, org):

	gff_name=strain_dict[name]  #need to open gff file 
	gff_name='gff/'+str(gff_name)

	gff_lib={} ##node:[[start,end,gene],[start,end,gene],....,] 
	gff_file=open(gff_name, 'r')

	for line in gff_file:


		if (line[0]!='#' and line[0]!='>' and line[0]!='A' and line[0]!='T' and line[0]!='G' and line[0]!='C'):
			line_list=line.strip().split('\t')


			info=str(line_list[8])

			if 'gene=' in info:
				temp=info.split('gene=')[1]
				gene=temp.split(';')[0]

			else:
				gene='hypothetical protein'

			temp=[]
			node=line_list[0]

				
			temp.append(line_list[3])
			temp.append(line_list[4])
			temp.append(gene) ##location of 'gene' in info tab too inconsistent, more over some don't even have it


		
			if node not in gff_lib:
				gff_lib[node]=[]
				gff_lib[node].append(temp)
			else:
				gff_lib[node].append(temp)

	gff_file.close()

	return gff_lib





##final function to mach between node LIST and node LIBRARY
def nodestogenes(nodes, gfflib):

	newdict={} #node_pos: gene to be used to translate  matrix


	
	for i in nodes:


		split_list=i.split('@')
		node=split_list[0]
		pos=int(split_list[1])
		
		repair=node+'_'+str(pos) ##original form for dictionary

		if node in gfflib:
			
			list_pos=gfflib[node]  ##list of [start, end, gene], [start, end, gene], 

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


			
			

			if gene=='': #because no matches were found at all
				#temp='intron'
				gene='intron'	


				"""
				for j in range(1000):
					pos2=pos+j
				
					pos3=pos-j
					

					if int(list_pos[k][0]) < pos2 and int(list_pos[k][1]) > pos2: # any hit 
						n+=1
						#m+=1 ##do wee need this???
						if n==1: ##first hit 
							temp=temp+'_ds'+str(j)+'_'+str(list_pos[k][2])
							break

						

					if int(list_pos[k][0]) < pos3 and int(list_pos[k][1]) > pos3: # any hit 
						n+=1

						if n==1: ##first hit 
							temp=temp+'_us'+str(j)+'_'+str(list_pos[k][2])
							break




				gene=temp
				"""

			
			newdict[repair]=gene

				

		else:
			#print(i, 'node not in dictionary')
			newdict[repair]='unknown node'

	#print('length of nodes list: ', len(nodes))
	#print('length of newdict: ', len(newdict))
	

	return(newdict)



def rewritematrix(patient):
	

	#######
	#get list of nodes 

	
	matrix_file='snp_matrix/'+ patient+'_merged_parsed.txt'


	mfile=open(matrix_file, 'r')
	nodelist=[]
	
	
	convertedmatrix={} ##{genes: [gene1,gene2,gene3,...], sample: [0,1,0,1,...]}
	convertedmatrix['genes']=[]

	first='T'
	for line in mfile:
		
		if first=='T':
			nodelist=line.strip().split('\t')
			for i in nodelist:

				gene=newdict[i]
				#print(i, gene)
				convertedmatrix['genes'].append(gene)

			first='F'

		elif first=='F': ##need to take info from matrix
			sampleline=line.strip().split('\t')
			
			sample=sampleline[0][48:]
			sample=sample[0:9]
			
			#print(sample)
			
			convertedmatrix[sample]=[]


			n =len(sampleline)
		
			for i in range(1,n):
				convertedmatrix[sample].append(sampleline[i])

	mfile.close()
	

	return convertedmatrix


def writefile(newmatrix):

	filename='snp_new_matrices/'+patient+'_merged_parsed.txt'
	newfile=open(filename, 'w')
	for i in newmatrix:
		n=len(newmatrix[i])
		
		newfile.write(i)
		newfile.write('\t')
		for k in range(n):
			newfile.write(newmatrix[i][k])
			newfile.write('\t')

		newfile.write('\n')


	newfile.close()



def countgenes(newmatrix):

	genes=newmatrix['genes']
	n = len(genes)
	m=len(newmatrix)

	temp={}
	introns={}

	for i in range(n):
		
		if 'intron_' not in genes[i]:
			if genes[i] not in temp:
				temp[genes[i]]=0

			for j in newmatrix:
				if j!='genes':
					temp[genes[i]] +=int(newmatrix[j][i])
	"""		

		if 'intron_' in genes[i]:

			intron_list=genes[i].strip().split('_')		
	
			gene2='_'.join(intron_list[2:])

			if gene2 in introns:
				introns[gene2].append(intron_list[1])
			else:
				introns[gene2]=[]
				introns[gene2].append(intron_list[1])

			
			'''
			if gene2 not in temp:
				temp[gene2]=0

			for j in newmatrix:
				if j!='genes':
					temp[gene2] +=int(newmatrix[j][i])
			'''
	"""
				
	return(temp)
	#return(temp, introns)


##############################################################################################################
#########################################EXECUTION FROM HERE ############################################
##############################################################################################################

if __name__ == '__main__':
    main()



n = len(matrices) #19
total=0

parallel={} #dict for counting # strains have a given mutation

for i in range(n):
	patient=matrices[i] ##matrix file name 
	
	#patient="MAB_20"
	name=patient ##+'_01'
	org='E_coli'

	gfflib=matchgff(name, org) #library of node:[[start,end,gene],[start,end,gene],....,] 
	nodes=matchnodes(patient) ##list of NODE@POS


	l=len(nodes) ##need to make sure node list is not 0 bc everything is 0 


	
	if l!=0:


		newdict=nodestogenes(nodes,gfflib)
		#print(newdict)
		##now that we have the dictionary, Write onto the header file of each matrix the new names 

		newmatrix=rewritematrix(patient)
		#print(newmatrix)

		temp=newmatrix['genes']
		m=len(temp)
		tempdict={} 

		for k in range(m):
			if temp[k] not in tempdict:
				tempdict[temp[k]]=1

		
		## optional : saves a record of which genes are mutated in parallel..
		for key in tempdict:
			if key not in parallel:
				parallel[key]=1
			else:
				parallel[key]+=1

		#writefile(newmatrix)


		t=countgenes(newmatrix)
	

		for item in t:
			print(patient, item, t[item], sep=',')
	
	

"""
for key in parallel:
	print(key, parallel[key], sep=',')

"""
