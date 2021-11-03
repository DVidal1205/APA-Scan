import os
import time


def SamtoText(input_path, bamfile, chromosomes):
	cwd = (os.getcwd()).split("/")[0:-1]
	cwd = '/'.join(cwd)
	output_dir = input_path+'/'+bamfile[:-4]+'/'
	os.makedirs(output_dir, exist_ok=True)
	cmd1 = cwd+"/samtools index "+input_path+"/"+bamfile		# make samtools index bamfile.bam.bai
	os.system(cmd1)
	print("Aligning bam files for",bamfile,"...")
	for chrom in chromosomes:
		tt = time.time()
		cmd2 = cwd+"/samtools view -b "+input_path+"/"+bamfile+" "+chrom+" -o "+output_dir+'/'+chrom+".bam"
		### Need to use pileup from SAMtools v0.1.8, not mpileup
		cmd3 = cwd+"/samtools pileup "+output_dir+'/'+chrom+".bam | cut -f 1,2,4 > "+output_dir+'/'+chrom+".txt" 
		command = cmd2+";"+cmd3
		os.system(command)
	return
	
