'''
------------------------------------------------------------------------------------------------------------------------
This workflow is used to run HybPiper on GenomeDK to investigate the paralogs in Ceroxyloids 
------------------------------------------------------------------------------------------------------------------------
This code is a variant of species_workflow.py by Oscar Wrisberg
Edited by Paola de Lima Ferreira 14/07/2022
------------------------------------------------------------------------------------------------------------------------
Eddited by Kristine Nørtoft Kristensen 
------------------------------------------------------------------------------------------------------------------------
'''
from os import O_SYNC, name  
from gwf import Workflow
import os.path   
import csv  

gwf = Workflow()

########################################################################################################################
################################################---- Hybpiper ----######################################################
########################################################################################################################
def hybpiper(sp, p1, p2, un, path_out, path_in, done):
    """Hybpiper with intronerate function.""" 
    inputs = [path_in + sp + p1, path_in + sp + p2, path_in + sp + un] # The files which the job will look for before it runs
    outputs = [path_out + sp, done] # The files which will have to be created in order for the job to be "completed"
    options = {'cores': 10, 'memory': "20g", 'walltime': "100:00:00", 'account':"bp_ceroxyloideae"} #Slurm commands 

    spec = """

    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate hybpiper

    cd /home/kris/bp_ceroxyloideae/hybpiper/
    hybpiper assemble --cpu 16 --readfiles {p1} {p2} --unpaired {un} --targetfile_dna /home/kris/palms_phylogeny/5.Ceroxyloids/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta  --prefix {sp} --bwa --run_intronerate
    
    cd /home/kris/bp_ceroxyloideae/done_hybpiper/
    touch {sp}
    
    """.format(sp=sp, p1 = path_in + sp + p1, p2 = path_in + sp + p2, un = path_in + sp + un , path_out = path_out, done = done)

    return (inputs, outputs, options, spec)

########################################################################################################################
###################################################---- Stats ----######################################################
#######################################################################################################################
def stats(path_in):
   """Gather statistics about the HybPiper run(s).""" 
   inputs = [path_in, "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero01", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero02", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero03", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero04", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero05", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero07", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero08", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero09", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero10", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero11", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero12", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero13","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero14", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero16", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero17", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero18", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero19", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero20", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero21", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero22", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero23", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero24", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero25", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero26", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero27", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero28", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero30", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero31", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero32", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero33", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero34", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero35", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero36", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero37", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero38", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero39", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero40", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero41", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero42", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero43", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero44", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero45", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero46", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero48", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero49", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero50", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero51", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero52", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero53", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero54", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero55", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero56", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero57", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero58", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero59", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero60", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero61", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero62", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero63", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero64", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero65", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero66", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero67", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero68", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero69", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero70", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero71","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero72","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero73","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero74","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero75","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero76","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero78","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero79","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero80","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero81","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero82", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero83","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero84", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero85", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero86","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero87", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero88", "/home/kris/bp_ceroxyloideae/done_hybpiper/Cero89","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero90","/home/kris/bp_ceroxyloideae/done_hybpiper/Cero91"]
   outputs = ["/home/kris/bp_ceroxyloideae/stats/"+"seq_lengths.tsv", "/home/kris/bp_ceroxyloideae/stats/"+"hybpiper_stats.tsv"]  # The files which will have to be created in order for the job to be "completed"
   options = {'cores': 10, 'memory': "20g", 'walltime': "100:00:00", 'account':"bp_ceroxyloideae"} #Slurm commands

   spec = """
   source /home/kris/miniconda3/etc/profile.d/conda.sh
   conda activate hybpiper
      
   cd /home/kris/bp_ceroxyloideae/hybpiper/
    
   hybpiper stats --targetfile_dna /home/kris/palms_phylogeny/5.Ceroxyloids/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta gene namelist.txt
     
   mv seq_lengths.tsv /home/kris/bp_ceroxyloideae/stats/
    
   mv hybpiper_stats.tsv /home/kris/bp_ceroxyloideae/stats/
      
   """.format(path_in = path_in)

   return (inputs, outputs, options, spec)


 

########################################################################################################################
#################################################---- Heatmap ----######################################################
########################################################################################################################

def heatmap(path_in):
   """Heatmap with Gene Recovery Information.""" 
   inputs = [path_in + "seq_lengths.tsv"] # The files which the job will look for before it runs
   outputs = ["/home/kris/bp_ceroxyloideae/heatmap/"+"recovery_heatmap.png"] # The files which will have to be created in order for the job to be "completed"
   options = {'cores': 10, 'memory': "20g", 'walltime': "24:00:00", 'account':"bp_ceroxyloideae"} #Slurm commands

   spec = """
   source /home/kris/miniconda3/etc/profile.d/conda.sh
   conda activate hybpiper
    
   cd /home/kris/bp_ceroxyloideae/stats/
    
   hybpiper recovery_heatmap seq_lengths.tsv
   
   mv recovery_heatmap* /home/kris/bp_ceroxyloideae/heatmap/
       
   """.format(path_in = path_in)

   return (inputs, outputs, options, spec)  
    
#Here you should stop and look at the graphic statistic for gene recovery 

########################################################################################################################
#############################################---- Paralogs 2 ----#######################################################
########################################################################################################################

def paralogs2(path_in):
   """Find Paralog genes and write them on the file called paralog.txt"""
   inputs = [path_in] #The files which the job will look for before it runs
   outputs = ["/home/kris/bp_ceroxyloideae/paralogs/"+"paralog_report.tsv", "/home/kris/bp_ceroxyloideae/paralogs/paralogs_all/", "/home/kris/bp_ceroxyloideae/paralogs/paralogs_no_chimeras/"]   # The files which will have to be created in order for the job to be "completed"
   options = {'cores': 2, 'memory': "10g", 'walltime': "8:00:00", 'account':"bp_ceroxyloideae"}

   spec = """
   source /home/kris/miniconda3/etc/profile.d/conda.sh 
   conda activate hybpiper

   cd {path_in}

   hybpiper paralog_retriever /home/kris/bp_ceroxyloideae/hybpiper/namelist.txt -t_dna /home/kris/palms_phylogeny/5.Ceroxyloids/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta 
   
   mv paralog_report.tsv /home/kris/bp_ceroxyloideae/paralogs/
   mv paralog_heatmap.png /home/kris/bp_ceroxyloideae/paralogs/
   mv paralogs_above_threshold_report.txt /home/kris/bp_ceroxyloideae/paralogs/
   mv paralogs_all/ /home/kris/bp_ceroxyloideae/paralogs/
   mv paralogs_no_chimeras/ /home/kris/bp_ceroxyloideae/paralogs/
 
   """.format(path_in = path_in)
   return(inputs, outputs, options, spec)



# ########################################################################################################################
# #############################################---- Coverage ----#########################################################
# ########################################################################################################################

#This script does the following:
#Gather all contigs from each sample in one fasta file: coverage/sample.fasta
#Map paired and unpaired reads to that fasta using BWA mem
#Deduplicate reads using Picard
#Calculate depth using samtools
#Mask/strip any bases with coverage <2
#Generate a new trimmed sample-level fasta: coverage/sample_trimmed.fasta


def coverage2(sp2, path_in, path_out, done,all_bam,all_sorted_bam, all_sorted_bam_bai, bam, cov,fasta,fasta_amb,fasta_ann,fasta_bwt,fasta_pac,fasta_sa,trimmed_fasta,up_bam,dir_in,dir_out,dir_wrk):
    """Calculating coverage of sequences."""
    inputs = [path_in+sp2]
    outputs = [path_out+sp2+all_bam, path_out+sp2+all_sorted_bam, path_out+sp2+all_sorted_bam_bai, path_out+sp2+bam,
    path_out+sp2+cov, path_out+sp2+fasta, path_out+sp2+fasta_amb, path_out+sp2+fasta_ann, path_out+sp2+fasta_bwt,
    path_out+sp2+fasta_pac, path_out+sp2+fasta_sa, path_out+sp2+trimmed_fasta, path_out+sp2+up_bam,done] #ALL the output files
    options = {'cores': 4, 'memory': "20g", 'walltime': "08:00:00", 'account':"bp_ceroxyloideae"}

    spec = """
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate coverage
    
    cd {path_in}

    python3.1 /home/kris/bp_ceroxyloideae/coverage/coverage.py {sample} {directory_in} {directory_out} {directory_wrk}

    touch {done}

    """.format(sample = sp2, done = done, path_in = path_in, directory_in = dir_in, directory_out = dir_out, directory_wrk = dir_wrk)

    return (inputs, outputs, options, spec)

    # ########################################################################################################################
# #############################################---- Retrieve 2 ----#######################################################
# ########################################################################################################################

#Think about doing blacklisting here? you could just remove species from the inputs here if you dont want them in the downstream analysis

def retrieve2(path_in):
    """Retrieve gene sequences from all the species and create an unaligned multifasta for each gene."""
    inputs = [path_in]
    outputs = ["/home/kris/bp_ceroxyloideae/coverage/Retrieve_all_done.txt"]
    options = {'cores': 10, 'memory': "20g", 'walltime': "06:00:00", 'account':"bp_ceroxyloideae"}

    spec = """

    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate coverage

    cd /home/kris/bp_ceroxyloideae/coverage/

    ls *trimmed.fasta > filelist.txt

    python3.1 sample2genes.py > outstats.csv

    touch /home/kris/bp_ceroxyloideae/coverage/Retrieve_all_done.txt

    """.format(path_in = path_in)

    return (inputs, outputs, options, spec)
    
###Here you should wait for the output. The output will comprise a file for each gene with the species sequence recovered.



##########################################################################################################################
###############################################---- MAFT ----#############################################################
##########################################################################################################################

# Here go to folder 6.Retrieve and ls -1
# Get the gene names and write them in genes = []
# We found 3467 genes for Ceroxyloids using the Arecoideae target file.

def mafft(genes, path_in, path_out, done):
    """Aligning all the sequences for each gene."""
    inputs = [path_in+genes]
    outputs = [done,path_out+genes+"_aligned.fasta"] 
    options = {'cores': 8, 'memory': "100g", 'walltime': "12:00:00", 'account':"bp_ceroxyloideae"}

    spec = """

    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate ame

    cd {path_in}

    mafft --auto --thread 64 {genes} > {path_out}{genes}_aligned.fasta
    
    touch {done}

    """.format(genes = genes, done = done, path_in = path_in, path_out=path_out)

    return (inputs, outputs, options, spec)

# ########################################################################################################################
# ###############################################---- TRIMAL ----#########################################################
# ########################################################################################################################

#Cleaning according trimmal
#Get raw alignments and trim them according to a gap threshold.
#Before you can run you need to mkdir 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95
#Yoy need a cutoff_trim.txt file with the cutoff listed, they should have the same value as the folders
#It is recommended to copy aligned.fasta files from mafft folder to a trimal folder

def gt_trimming(path_in,done,genes):
    """ Use trimmal for trimming all alignments for each of the GT values specified"""
    inputs = [path_in+genes+"_aligned.fasta"]
    outputs = [done]
    options = {'cores': 1, 'memory': "20g", 'walltime': "12:00:00", 'account':"bp_ceroxyloideae"}

    spec="""
    #Activating enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate trimal

    #Go to alignments folder
    cd /home/kris/bp_ceroxyloideae/trimal/

    #Running gaptrimming.sh
    
    trimal -in {genes}_aligned.fasta -out 0.1/{genes}_aligned.fasta -htmlout 0.1/{genes}_aligned.fasta.html -gt 0.1
    trimal -in {genes}_aligned.fasta -out 0.15/{genes}_aligned.fasta -htmlout 0.15/{genes}_aligned.fasta.html -gt 0.15
    trimal -in {genes}_aligned.fasta -out 0.2/{genes}_aligned.fasta -htmlout 0.2/{genes}_aligned.fasta.html -gt 0.2
    trimal -in {genes}_aligned.fasta -out 0.25/{genes}_aligned.fasta -htmlout 0.25/{genes}_aligned.fasta.html -gt 0.25
    trimal -in {genes}_aligned.fasta -out 0.3/{genes}_aligned.fasta -htmlout 0.3/{genes}_aligned.fasta.html -gt 0.3
    trimal -in {genes}_aligned.fasta -out 0.35/{genes}_aligned.fasta -htmlout 0.35/{genes}_aligned.fasta.html -gt 0.35
    trimal -in {genes}_aligned.fasta -out 0.4/{genes}_aligned.fasta -htmlout 0.4/{genes}_aligned.fasta.html -gt 0.4
    trimal -in {genes}_aligned.fasta -out 0.45/{genes}_aligned.fasta -htmlout 0.45/{genes}_aligned.fasta.html -gt 0.45
    trimal -in {genes}_aligned.fasta -out 0.5/{genes}_aligned.fasta -htmlout 0.5/{genes}_aligned.fasta.html -gt 0.5
    trimal -in {genes}_aligned.fasta -out 0.55/{genes}_aligned.fasta -htmlout 0.55/{genes}_aligned.fasta.html -gt 0.55
    trimal -in {genes}_aligned.fasta -out 0.6/{genes}_aligned.fasta -htmlout 0.6/{genes}_aligned.fasta.html -gt 0.6
    trimal -in {genes}_aligned.fasta -out 0.65/{genes}_aligned.fasta -htmlout 0.65/{genes}_aligned.fasta.html -gt 0.65
    trimal -in {genes}_aligned.fasta -out 0.7/{genes}_aligned.fasta -htmlout 0.7/{genes}_aligned.fasta.html -gt 0.7
    trimal -in {genes}_aligned.fasta -out 0.75/{genes}_aligned.fasta -htmlout 0.75/{genes}_aligned.fasta.html -gt 0.75
    trimal -in {genes}_aligned.fasta -out 0.8/{genes}_aligned.fasta -htmlout 0.8/{genes}_aligned.fasta.html -gt 0.8
    trimal -in {genes}_aligned.fasta -out 0.85/{genes}_aligned.fasta -htmlout 0.85/{genes}_aligned.fasta.html -gt 0.85
    trimal -in {genes}_aligned.fasta -out 0.9/{genes}_aligned.fasta -htmlout 0.9/{genes}_aligned.fasta.html -gt 0.9
    trimal -in {genes}_aligned.fasta -out 0.95/{genes}_aligned.fasta -htmlout 0.95/{genes}_aligned.fasta.html -gt 0.95
    
    cat 0.*/{genes}_aligned.fasta.html > /home/kris/bp_ceroxyloideae/done_trimal/{genes}_aligned.fasta.html_done
    
    touch {done}

    """.format(path_in=path_in, done=done, genes=genes)

    return(inputs, outputs, options, spec)  




# ########################################################################################################################
# #################################################---- AMAS ----#########################################################
# ########################################################################################################################

######## Calculating amas summary
# For raw alignments
def amas_raw(path_in):
    """Creating summary files for all the trimmed alignments for each raw alignment"""
    inputs = [path_in]
    outputs = [path_in+"summary_0.txt",path_in+"summary_0.1.txt",path_in+"summary_0.15.txt",path_in+"summary_0.2.txt",path_in+"summary_0.25.txt",path_in+"summary_0.3.txt", path_in+"summary_0.35.txt",path_in+"summary_0.4.txt",path_in+"summary_0.45.txt",path_in+"summary_0.5.txt",path_in+"summary_0.55.txt",path_in+"summary_0.6.txt",path_in+"summary_0.65.txt", path_in+"summary_0.7.txt",path_in+"summary_0.75.txt",path_in+"summary_0.8.txt",path_in+"summary_0.85.txt",path_in+"summary_0.9.txt",path_in+"summary_0.95.txt"]
    options = {'cores': 1, 'memory': "10g", 'walltime': "12:00:00", 'account':"bp_ceroxyloideae"}

    spec="""

    #Activating enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate amas
    
    #Calculating amas summary
    bash /home/kris/bp_ceroxyloideae/scripts/amas_raw.sh
    bash /home/kris/bp_ceroxyloideae/scripts/amas_gt1.sh
    
    """.format(path_in = path_in)

    return(inputs, outputs, options, spec)


#Hereafter you need to do some manual work and ad headlines statistics from one of the gene samples to the summary files.


    
# ########################################################################################################################
# #############################################---- Optrimal ----#########################################################
# ########################################################################################################################

#Getting the best alignment for each gene 
def optrim(path_in):
    """Select the best alignments according to the gt value"""
    inputs = [path_in+"summary_0.txt",path_in+"summary_0.1.txt",path_in+"summary_0.15.txt",path_in+"summary_0.2.txt",path_in+"summary_0.25.txt",path_in+"summary_0.3.txt",
    path_in+"summary_0.35.txt",path_in+"summary_0.4.txt",path_in+"summary_0.45.txt",path_in+"summary_0.5.txt",path_in+"summary_0.55.txt",path_in+"summary_0.6.txt",path_in+"summary_0.65.txt",
    path_in+"summary_0.7.txt",path_in+"summary_0.75.txt",path_in+"summary_0.8.txt",path_in+"summary_0.85.txt",path_in+"summary_0.9.txt",path_in+"summary_0.95.txt"]
    outputs = ["/home/kris/bp_ceroxyloideae/best_alignments/optimal_final_results/"]
    options = {'cores': 10, 'memory': "20g", 'walltime': "08:00:00", 'account':"bp_ceroxyloideae"}

    spec="""

    #Going to folder with trimmed files
    cd {path_in}

    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate R
    
    Rscript --vanilla /home/kris/bp_ceroxyloideae/scripts/optrimal.R

    mv dldp_* optrim_output/
    
    cd /home/kris/bp_ceroxyloideae/best_alignments/optimal_final_results 
    
    cp * /home/kris/bp_ceroxyloideae/best_alignments

    """.format(path_in=path_in)

    return(inputs, outputs, options, spec)



##########################################################################################################################
##############################################---- CIALIGN ----###########################################################
##########################################################################################################################

def cialign1(genes, path_in, path_out):
    """Cleaning alignments using cialign default."""
    inputs = [path_in + genes + "_aligned.fasta"]
    outputs = [path_out+genes+"_cialign.fasta_cleaned.fasta"]
    options = {'cores': 8, 'memory': "100g", 'walltime': "12:00:00", 'account':"bp_ceroxyloideae"}

    spec = """
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate CIAlign

    cd /home/kris/bp_ceroxyloideae/best_alignments/optimal_final_results/

    CIAlign --infile {genes}_aligned.fasta --all --outfile_stem {path_out}{genes}_cialign.fasta

    """.format(genes = genes, path_in = path_in, path_out=path_out)

    return (inputs, outputs, options, spec)


# ########################################################################################################################
# ###############################################---- TAPER ----##########################################################
# ########################################################################################################################

def taper(path_in, genes, path_out):
    """Using TAPER AFTER CIAlign to remove errors in small species-specific stretches of the multiple sequence alignments"""
    inputs = [path_in+genes+"_cialign.fasta_cleaned.fasta"]
    outputs = ["/home/kris/bp_ceroxyloideae/taper/"+genes+"_output_tapper.fasta"]
    options = {'cores': 1, 'memory': "40g", 'walltime': "02:00:00", 'account':"bp_ceroxyloideae"}

    spec = """
     
    cd /home/kris/bp_ceroxyloideae/cialign/
        
    #Activate the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate julia
        
    julia /home/kris/miniconda3/envs/julia/correction_multi.jl {genes}_cialign.fasta_cleaned.fasta > {genes}_output_tapper.fasta 
    

    mv *_output_tapper.fasta /home/kris/bp_ceroxyloideae/taper/
        
    """.format(path_in = path_in, genes = genes, path_out = path_out)

    return (inputs, outputs, options, spec)

# ########################################################################################################################
# ##############################################---- IQTREE ----##########################################################
# ########################################################################################################################

def iqtree(path_in, genes):
    """Using IQTREE to construct a phylogenetic hypotheses for each gene"""
    inputs = [path_in+genes+"_output_tapper.fasta"]
    outputs = ["/home/kris/bp_ceroxyloideae/IQtree/"+genes+"_output_tapper.fasta.treefile"]
    options = {'cores': 2, 'memory': "40g", 'walltime': "12:00:00", 'account':"bp_ceroxyloideae"}

    spec = """
     
    cd /home/kris/bp_ceroxyloideae/taper/
        
    #Activate the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate iqtree
        
    iqtree2 -s {genes}_output_tapper.fasta -T AUTO -m MFP -B 1000 
    
    mv *treefile /home/kris/bp_ceroxyloideae/IQtree
    mv *_output_tapper.fasta.model.gz /home/kris/bp_ceroxyloideae/IQtree
    mv *output_tapper.fasta.contree /home/kris/bp_ceroxyloideae/IQtree
    mv *output_tapper.fasta.bionj /home/kris/bp_ceroxyloideae/IQtree
    mv *output_tapper.fasta.ckp.gz /home/kris/bp_ceroxyloideae/IQtree
    mv *_output_tapper.fasta.iqtree /home/kris/bp_ceroxyloideae/IQtree
    mv *_output_tapper.fasta.log /home/kris/bp_ceroxyloideae/IQtree
    mv *tapper.fasta.mldist /home/kris/bp_ceroxyloideae/IQtree
    mv *tapper.fasta.splits.nex /home/kris/bp_ceroxyloideae/IQtree
    mv *output_tapper.fasta.uniqueseq.phy /home/kris/bp_ceroxyloideae/IQtree
    
    """.format(path_in = path_in, genes = genes)

    return (inputs, outputs, options, spec) 
 
#Here you should stop and go to folder /IQtree and do cat *treefile > gene_trees.nex


# ########################################################################################################################
# #####################################---- Astral Tree Search ----#####################################################
# ########################################################################################################################


#Showing Posterior Probabilities
def astral_2(path_in, gene_tree_file, output, genes):
    """Using Astral to construct a species tree based on the genetrees"""
    inputs = [path_in]
    outputs = [path_in + output]
    options = {'cores': 10, 'memory': "100g", 'walltime': "12:00:00", 'account':"bp_ceroxyloideae"}

    spec = """
    #source /home/paola/miniconda2/etc/profile.d/conda.sh
     
    cd /home/kris/Astral/
  
    java -D"java.library.path=lib/" -jar astral.5.7.8.jar -i /home/kris/bp_ceroxyloideae/astral/gene_trees.nex -o astral_tree_probabilities.tre 2> log_posterior_probability.out
    
    mv astral_tree_probabilities.tre /home/kris/bp_ceroxyloideae/astral/
    mv log_posterior_probability.out /home/kris/bp_ceroxyloideae/astral/
    """.format(path_in = path_in, gene_tree_file = gene_tree_file, output=output, genes = genes)
    
    return (inputs, outputs, options, spec)
        
#Showing number of genes supporting clades   
def astral_gene_supporting(path_in, gene_tree_file, output, genes):
    """Using Astral to construct a species tree based on the genetrees"""
    inputs = [path_in+"gene_trees.nex"]
    outputs = [path_in + output]
    options = {'cores': 10, 'memory': "40g", 'walltime': "24:00:00", 'account':"bp_ceroxyloideae"}

    spec = """
    #source /home/paola/miniconda2/etc/profile.d/conda.sh
                     
    cd /home/kris/Astral/ 
  
    java -D"java.library.path=lib/" -jar astral.5.7.8.jar -i /home/kris/bp_ceroxyloideae/astral/gene_trees.nex -o astral_gene_probability.tre -t 2 2> log_gene_effective.out
      
    mv astral_gene_probability.tre /home/kris/bp_ceroxyloideae/astral/
    mv log_gene_effective.out /home/kris/bp_ceroxyloideae/astral/
    
    """.format(path_in = path_in, gene_tree_file = gene_tree_file, output=output, genes = genes)

    return (inputs, outputs, options, spec)



########################################################################################################################
#####################################################---- RUN ----#######################################################
########################################################################################################################

sp = ["Cero01",  "Cero02",  "Cero03",  "Cero04",  "Cero05",  "Cero07",  "Cero08",  "Cero09",  "Cero10",  "Cero11",  "Cero12",  "Cero13",  "Cero14",  "Cero16",  "Cero17",  "Cero18",  "Cero19",  "Cero20",  "Cero21",  "Cero22",  "Cero23",  "Cero24",  "Cero25",  "Cero26",  "Cero27",  "Cero28",  "Cero30",  "Cero31",  "Cero32",  "Cero33",  "Cero34",  "Cero35",  "Cero36",  "Cero37",  "Cero38",  "Cero39",  "Cero40",  "Cero41",  "Cero42",  "Cero43",  "Cero44",  "Cero45",  "Cero46",  "Cero48",  "Cero49",  "Cero50",  "Cero51",  "Cero52",  "Cero53",  "Cero54",  "Cero55",  "Cero56",  "Cero57",  "Cero58",  "Cero59",  "Cero60",  "Cero61",  "Cero62",  "Cero63",  "Cero64",  "Cero65",  "Cero66",  "Cero67",  "Cero68",  "Cero69",  "Cero70",  "Cero71", "Cero72", "Cero73", "Cero74", "Cero75", "Cero76", "Cero77", "Cero78", "Cero79", "Cero80", "Cero81", "Cero82", "Cero83", "Cero84", "Cero85", "Cero86", "Cero87", "Cero88", "Cero89", "Cero90", "Cero91"]

genes=["EGU105032175.FNA","EGU105032229.FNA","EGU105032337.FNA","EGU105032379.FNA","EGU105033063.FNA","EGU105033626.FNA","EGU105034121.FNA","EGU105034616.FNA","EGU105034893.FNA",
"EGU105034993.FNA","EGU105035046.FNA","EGU105035196.FNA","EGU105035203.FNA","EGU105035462.FNA","EGU105035555.FNA","EGU105035989.FNA","EGU105036031.FNA","EGU105036385.FNA",
"EGU105036774.FNA","EGU105037749.FNA","EGU105037800.FNA","EGU105037890.FNA","EGU105037902.FNA","EGU105037930.FNA","EGU105037938.FNA","EGU105038008.FNA",
"EGU105038036.FNA","EGU105038098.FNA","EGU105038099.FNA","EGU105038100.FNA","EGU105038110.FNA","EGU105038114.FNA","EGU105038118.FNA","EGU105038123.FNA",
"EGU105038179.FNA","EGU105038201.FNA","EGU105038228.FNA","EGU105038234.FNA","EGU105038245.FNA","EGU105038252.FNA","EGU105038310.FNA","EGU105038382.FNA",
"EGU105038400.FNA","EGU105038419.FNA","EGU105038431.FNA","EGU105038499.FNA","EGU105038513.FNA","EGU105038571.FNA","EGU105038580.FNA","EGU105038603.FNA",
"EGU105038631.FNA","EGU105038680.FNA","EGU105038693.FNA","EGU105038720.FNA","EGU105038747.FNA","EGU105038794.FNA","EGU105038832.FNA","EGU105038882.FNA",
"EGU105038986.FNA","EGU105038988.FNA","EGU105039013.FNA","EGU105039062.FNA","EGU105039067.FNA","EGU105039082.FNA","EGU105039099.FNA","EGU105039101.FNA",
"EGU105039107.FNA","EGU105039121.FNA","EGU105039164.FNA","EGU105039178.FNA","EGU105039221.FNA","EGU105039236.FNA","EGU105039255.FNA",
"EGU105039282.FNA","EGU105039298.FNA","EGU105039313.FNA","EGU105039403.FNA","EGU105039431.FNA","EGU105039449.FNA","EGU105039460.FNA","EGU105039480.FNA",
"EGU105039494.FNA","EGU105039501.FNA","EGU105039512.FNA","EGU105039542.FNA","EGU105039587.FNA","EGU105039595.FNA","EGU105039609.FNA","EGU105039660.FNA",
"EGU105039685.FNA","EGU105039690.FNA","EGU105039699.FNA","EGU105039763.FNA","EGU105039783.FNA","EGU105039809.FNA","EGU105039822.FNA","EGU105039925.FNA",
"EGU105039947.FNA","EGU105039957.FNA","EGU105039962.FNA","EGU105040073.FNA","EGU105040088.FNA","EGU105040099.FNA","EGU105040114.FNA","EGU105040115.FNA","EGU105040125.FNA",
"EGU105040139.FNA","EGU105040185.FNA","EGU105040186.FNA","EGU105040189.FNA","EGU105040206.FNA","EGU105040207.FNA","EGU105040242.FNA","EGU105040281.FNA",
"EGU105040302.FNA","EGU105040308.FNA","EGU105040359.FNA","EGU105040368.FNA","EGU105040426.FNA","EGU105040452.FNA","EGU105040462.FNA","EGU105040530.FNA",
"EGU105040583.FNA","EGU105040667.FNA","EGU105040675.FNA","EGU105040684.FNA","EGU105040690.FNA","EGU105040700.FNA","EGU105040756.FNA","EGU105040758.FNA",
"EGU105040813.FNA","EGU105040837.FNA","EGU105040842.FNA","EGU105040850.FNA","EGU105040851.FNA","EGU105040863.FNA","EGU105040887.FNA","EGU105040914.FNA",
"EGU105040918.FNA","EGU105040922.FNA","EGU105040957.FNA","EGU105040970.FNA","EGU105041055.FNA","EGU105041100.FNA","EGU105041117.FNA","EGU105041125.FNA",
"EGU105041127.FNA","EGU105041133.FNA","EGU105041179.FNA","EGU105041182.FNA","EGU105041189.FNA","EGU105041217.FNA","EGU105041246.FNA","EGU105041283.FNA",
"EGU105041337.FNA","EGU105041353.FNA","EGU105041650.FNA","EGU105041657.FNA","EGU105041665.FNA","EGU105041680.FNA","EGU105041687.FNA","EGU105041710.FNA",
"EGU105041807.FNA","EGU105041816.FNA","EGU105041872.FNA","EGU105041902.FNA","EGU105041903.FNA","EGU105041929.FNA","EGU105041933.FNA","EGU105041982.FNA",
"EGU105042090.FNA","EGU105042113.FNA","EGU105042128.FNA","EGU105042147.FNA","EGU105042168.FNA","EGU105042205.FNA","EGU105042290.FNA","EGU105042307.FNA",
"EGU105042323.FNA","EGU105042329.FNA","EGU105042368.FNA","EGU105042422.FNA","EGU105042525.FNA","EGU105042558.FNA","EGU105042560.FNA","EGU105042584.FNA",
"EGU105042633.FNA","EGU105042644.FNA","EGU105042651.FNA","EGU105042664.FNA","EGU105042722.FNA","EGU105042781.FNA","EGU105042808.FNA","EGU105042820.FNA",
"EGU105042873.FNA","EGU105042965.FNA","EGU105043011.FNA","EGU105043037.FNA","EGU105043042.FNA","EGU105043061.FNA","EGU105043069.FNA",
"EGU105043119.FNA","EGU105043155.FNA","EGU105043160.FNA","EGU105043164.FNA","EGU105043193.FNA","EGU105043320.FNA","EGU105043338.FNA","EGU105043374.FNA",
"EGU105043419.FNA","EGU105043430.FNA","EGU105043469.FNA","EGU105043485.FNA","EGU105043499.FNA","EGU105043601.FNA","EGU105043633.FNA","EGU105043666.FNA",
"EGU105043685.FNA","EGU105043686.FNA","EGU105043730.FNA","EGU105043786.FNA","EGU105043816.FNA","EGU105043827.FNA","EGU105043926.FNA","EGU105043975.FNA",
"EGU105044063.FNA","EGU105044120.FNA","EGU105044133.FNA","EGU105044174.FNA","EGU105044182.FNA","EGU105044183.FNA","EGU105044203.FNA",
"EGU105044252.FNA","EGU105044281.FNA","EGU105044307.FNA","EGU105044309.FNA","EGU105044350.FNA","EGU105044378.FNA","EGU105044400.FNA","EGU105044407.FNA",
"EGU105044445.FNA","EGU105044446.FNA","EGU105044481.FNA","EGU105044588.FNA","EGU105044613.FNA","EGU105044614.FNA","EGU105044668.FNA",
"EGU105044676.FNA","EGU105044710.FNA","EGU105044758.FNA","EGU105044844.FNA","EGU105044846.FNA","EGU105044854.FNA","EGU105044885.FNA","EGU105044893.FNA",
"EGU105044896.FNA","EGU105044978.FNA","EGU105044982.FNA","EGU105044983.FNA","EGU105044984.FNA","EGU105045005.FNA","EGU105045043.FNA","EGU105045070.FNA",
"EGU105045078.FNA","EGU105045094.FNA","EGU105045099.FNA","EGU105045102.FNA","EGU105045137.FNA","EGU105045148.FNA","EGU105045232.FNA","EGU105045248.FNA",
"EGU105045254.FNA","EGU105045282.FNA","EGU105045310.FNA","EGU105045358.FNA","EGU105045367.FNA","EGU105045424.FNA","EGU105045464.FNA","EGU105045467.FNA",
"EGU105045489.FNA","EGU105045507.FNA","EGU105045509.FNA","EGU105045514.FNA","EGU105045520.FNA","EGU105045529.FNA","EGU105045544.FNA","EGU105045640.FNA",
"EGU105045658.FNA","EGU105045703.FNA","EGU105045726.FNA","EGU105045732.FNA","EGU105045760.FNA","EGU105045782.FNA","EGU105045788.FNA","EGU105045820.FNA",
"EGU105045827.FNA","EGU105045828.FNA","EGU105045835.FNA","EGU105045898.FNA","EGU105045932.FNA","EGU105045946.FNA","EGU105046030.FNA","EGU105046050.FNA",
"EGU105046056.FNA","EGU105046099.FNA","EGU105046103.FNA","EGU105046147.FNA","EGU105046168.FNA","EGU105046245.FNA","EGU105046297.FNA","EGU105046360.FNA",
"EGU105046387.FNA","EGU105046393.FNA","EGU105046401.FNA","EGU105046449.FNA","EGU105046454.FNA","EGU105046456.FNA","EGU105046503.FNA","EGU105046518.FNA",
"EGU105046530.FNA","EGU105046549.FNA","EGU105046559.FNA","EGU105046562.FNA","EGU105046574.FNA","EGU105046630.FNA","EGU105046632.FNA","EGU105046696.FNA",
"EGU105046735.FNA","EGU105046766.FNA","EGU105046786.FNA","EGU105046827.FNA","EGU105046875.FNA","EGU105046918.FNA","EGU105047024.FNA","EGU105047029.FNA",
"EGU105047253.FNA","EGU105047288.FNA","EGU105047293.FNA","EGU105047342.FNA","EGU105047357.FNA","EGU105047362.FNA","EGU105047379.FNA","EGU105047385.FNA",
"EGU105047395.FNA","EGU105047433.FNA","EGU105047434.FNA","EGU105047446.FNA","EGU105047519.FNA","EGU105047533.FNA","EGU105047546.FNA","EGU105047553.FNA",
"EGU105047578.FNA","EGU105047585.FNA","EGU105047597.FNA","EGU105047621.FNA","EGU105047644.FNA","EGU105047662.FNA","EGU105047689.FNA","EGU105047751.FNA",
"EGU105047777.FNA","EGU105047790.FNA","EGU105047907.FNA","EGU105047916.FNA","EGU105047922.FNA","EGU105047940.FNA","EGU105047945.FNA","EGU105047970.FNA",
"EGU105048009.FNA","EGU105048015.FNA","EGU105048028.FNA","EGU105048054.FNA","EGU105048056.FNA","EGU105048129.FNA","EGU105048130.FNA","EGU105048137.FNA",
"EGU105048159.FNA","EGU105048182.FNA","EGU105048199.FNA","EGU105048300.FNA","EGU105048357.FNA","EGU105048410.FNA","EGU105048474.FNA","EGU105048476.FNA",
"EGU105048479.FNA","EGU105048484.FNA","EGU105048486.FNA","EGU105048493.FNA","EGU105048527.FNA","EGU105048541.FNA","EGU105048581.FNA","EGU105048612.FNA",
"EGU105048694.FNA","EGU105048725.FNA","EGU105048751.FNA","EGU105048796.FNA","EGU105048839.FNA","EGU105048867.FNA","EGU105048886.FNA","EGU105048898.FNA",
"EGU105048909.FNA","EGU105048915.FNA","EGU105048926.FNA","EGU105048961.FNA","EGU105048968.FNA","EGU105049007.FNA","EGU105049016.FNA","EGU105049020.FNA",
"EGU105049025.FNA","EGU105049052.FNA","EGU105049097.FNA","EGU105049274.FNA","EGU105049312.FNA","EGU105049318.FNA","EGU105049360.FNA","EGU105049426.FNA",
"EGU105049539.FNA","EGU105049543.FNA","EGU105049583.FNA","EGU105049690.FNA","EGU105049729.FNA","EGU105049737.FNA","EGU105049761.FNA","EGU105049827.FNA",
"EGU105049882.FNA","EGU105049902.FNA","EGU105049903.FNA","EGU105049934.FNA","EGU105049947.FNA","EGU105050012.FNA","EGU105050023.FNA","EGU105050036.FNA",
"EGU105050058.FNA","EGU105050114.FNA","EGU105050126.FNA","EGU105050202.FNA","EGU105050207.FNA","EGU105050328.FNA","EGU105050344.FNA","EGU105050362.FNA",
"EGU105050366.FNA","EGU105050383.FNA","EGU105050387.FNA","EGU105050404.FNA","EGU105050432.FNA","EGU105050450.FNA","EGU105050521.FNA","EGU105050532.FNA",
"EGU105050644.FNA","EGU105050670.FNA","EGU105050680.FNA","EGU105050681.FNA","EGU105050682.FNA","EGU105050831.FNA","EGU105050841.FNA","EGU105050853.FNA",
"EGU105050854.FNA","EGU105050961.FNA","EGU105050970.FNA","EGU105050972.FNA","EGU105051087.FNA","EGU105051146.FNA","EGU105051156.FNA","EGU105051188.FNA",
"EGU105051345.FNA","EGU105051362.FNA","EGU105051366.FNA","EGU105051373.FNA","EGU105051391.FNA","EGU105051395.FNA","EGU105051403.FNA","EGU105051481.FNA",
"EGU105051499.FNA","EGU105051503.FNA","EGU105051560.FNA","EGU105051564.FNA","EGU105051582.FNA","EGU105051614.FNA","EGU105051677.FNA","EGU105051704.FNA",
"EGU105051726.FNA","EGU105051740.FNA","EGU105051748.FNA","EGU105051764.FNA","EGU105051795.FNA","EGU105051802.FNA","EGU105051821.FNA","EGU105051823.FNA",
"EGU105051832.FNA","EGU105051847.FNA","EGU105051857.FNA","EGU105051860.FNA","EGU105051870.FNA","EGU105051891.FNA","EGU105051924.FNA","EGU105051953.FNA",
"EGU105051985.FNA","EGU105052035.FNA","EGU105052070.FNA","EGU105052170.FNA","EGU105052178.FNA","EGU105052304.FNA","EGU105052307.FNA","EGU105052346.FNA",
"EGU105052351.FNA","EGU105052386.FNA","EGU105052389.FNA","EGU105052394.FNA","EGU105052428.FNA","EGU105052446.FNA","EGU105052476.FNA","EGU105052483.FNA",
"EGU105052492.FNA","EGU105052495.FNA","EGU105052527.FNA","EGU105052529.FNA","EGU105052538.FNA","EGU105052552.FNA","EGU105052573.FNA","EGU105052580.FNA",
"EGU105052623.FNA","EGU105052650.FNA","EGU105052694.FNA","EGU105052704.FNA","EGU105052739.FNA","EGU105052743.FNA","EGU105052750.FNA","EGU105052771.FNA",
"EGU105052804.FNA","EGU105052818.FNA","EGU105052849.FNA","EGU105052855.FNA","EGU105052865.FNA","EGU105052888.FNA","EGU105052944.FNA","EGU105052947.FNA",
"EGU105052956.FNA","EGU105053006.FNA","EGU105053055.FNA","EGU105053059.FNA","EGU105053079.FNA","EGU105053105.FNA","EGU105053124.FNA","EGU105053130.FNA",
"EGU105053136.FNA","EGU105053172.FNA","EGU105053204.FNA","EGU105053227.FNA","EGU105053263.FNA","EGU105053403.FNA","EGU105053422.FNA","EGU105053426.FNA",
"EGU105053457.FNA","EGU105053465.FNA","EGU105053468.FNA","EGU105053482.FNA","EGU105053549.FNA","EGU105053642.FNA","EGU105053654.FNA","EGU105053735.FNA",
"EGU105053747.FNA","EGU105053770.FNA","EGU105053835.FNA","EGU105053848.FNA","EGU105053866.FNA","EGU105053889.FNA","EGU105053901.FNA","EGU105053932.FNA",
"EGU105053961.FNA","EGU105053969.FNA","EGU105053974.FNA","EGU105053980.FNA","EGU105054002.FNA","EGU105054124.FNA","EGU105054130.FNA","EGU105054153.FNA",
"EGU105054204.FNA","EGU105054280.FNA","EGU105054293.FNA","EGU105054405.FNA","EGU105054435.FNA","EGU105054440.FNA","EGU105054455.FNA","EGU105054457.FNA","EGU105054469.FNA",
"EGU105054478.FNA","EGU105054486.FNA","EGU105054498.FNA","EGU105054529.FNA","EGU105054534.FNA","EGU105054595.FNA","EGU105054649.FNA","EGU105054653.FNA",
"EGU105054668.FNA","EGU105054723.FNA","EGU105054765.FNA","EGU105054786.FNA","EGU105054827.FNA","EGU105054845.FNA","EGU105054864.FNA","EGU105054891.FNA",
"EGU105054896.FNA","EGU105054898.FNA","EGU105054924.FNA","EGU105054930.FNA","EGU105054936.FNA","EGU105054948.FNA","EGU105054972.FNA","EGU105055008.FNA",
"EGU105055015.FNA","EGU105055023.FNA","EGU105055024.FNA","EGU105055030.FNA","EGU105055047.FNA","EGU105055052.FNA","EGU105055065.FNA","EGU105055072.FNA",
"EGU105055075.FNA","EGU105055077.FNA","EGU105055090.FNA","EGU105055093.FNA","EGU105055098.FNA","EGU105055114.FNA","EGU105055115.FNA","EGU105055130.FNA",
"EGU105055144.FNA","EGU105055157.FNA","EGU105055201.FNA","EGU105055283.FNA","EGU105055433.FNA","EGU105055438.FNA","EGU105055490.FNA","EGU105055499.FNA",
"EGU105055507.FNA","EGU105055550.FNA","EGU105055569.FNA","EGU105055621.FNA","EGU105055634.FNA","EGU105055664.FNA","EGU105055709.FNA","EGU105055755.FNA",
"EGU105055761.FNA","EGU105055771.FNA","EGU105055800.FNA","EGU105055862.FNA","EGU105055873.FNA","EGU105055883.FNA","EGU105055889.FNA","EGU105055908.FNA",
"EGU105055912.FNA","EGU105055913.FNA","EGU105056032.FNA","EGU105056091.FNA","EGU105056151.FNA","EGU105056269.FNA","EGU105056287.FNA","EGU105056289.FNA",
"EGU105056313.FNA","EGU105056323.FNA","EGU105056365.FNA","EGU105056382.FNA","EGU105056393.FNA","EGU105056460.FNA","EGU105056468.FNA","EGU105056469.FNA",
"EGU105056496.FNA","EGU105056530.FNA","EGU105056534.FNA","EGU105056539.FNA","EGU105056654.FNA","EGU105056662.FNA","EGU105056684.FNA","EGU105056688.FNA",
"EGU105056714.FNA","EGU105056726.FNA","EGU105056817.FNA","EGU105056848.FNA","EGU105056881.FNA","EGU105056943.FNA","EGU105056960.FNA","EGU105056998.FNA",
"EGU105057013.FNA","EGU105057015.FNA","EGU105057019.FNA","EGU105057074.FNA","EGU105057090.FNA","EGU105057110.FNA","EGU105057130.FNA","EGU105057194.FNA",
"EGU105057235.FNA","EGU105057256.FNA","EGU105057335.FNA","EGU105057357.FNA","EGU105057553.FNA","EGU105057579.FNA","EGU105057634.FNA","EGU105057666.FNA",
"EGU105057669.FNA","EGU105057721.FNA","EGU105057742.FNA","EGU105057795.FNA","EGU105057841.FNA","EGU105057912.FNA","EGU105057919.FNA","EGU105057941.FNA",
"EGU105058078.FNA","EGU105058081.FNA","EGU105058083.FNA","EGU105058094.FNA","EGU105058107.FNA","EGU105058131.FNA","EGU105058170.FNA","EGU105058175.FNA",
"EGU105058180.FNA","EGU105058202.FNA","EGU105058237.FNA","EGU105058241.FNA","EGU105058245.FNA","EGU105058326.FNA","EGU105058366.FNA","EGU105058377.FNA",
"EGU105058418.FNA","EGU105058469.FNA","EGU105058499.FNA","EGU105058547.FNA","EGU105058556.FNA","EGU105058567.FNA","EGU105058576.FNA","EGU105058582.FNA",
"EGU105058592.FNA","EGU105058598.FNA","EGU105058614.FNA","EGU105058633.FNA","EGU105058682.FNA","EGU105058683.FNA","EGU105058687.FNA","EGU105058702.FNA",
"EGU105058723.FNA","EGU105058731.FNA","EGU105058781.FNA","EGU105058798.FNA","EGU105058802.FNA","EGU105058808.FNA","EGU105058863.FNA","EGU105058889.FNA",
"EGU105058890.FNA","EGU105058894.FNA","EGU105058904.FNA","EGU105058918.FNA","EGU105058989.FNA","EGU105058990.FNA","EGU105059003.FNA","EGU105059008.FNA",
"EGU105059023.FNA","EGU105059035.FNA","EGU105059042.FNA","EGU105059054.FNA","EGU105059108.FNA","EGU105059112.FNA","EGU105059113.FNA","EGU105059126.FNA",
"EGU105059131.FNA","EGU105059138.FNA","EGU105059176.FNA","EGU105059186.FNA","EGU105059193.FNA","EGU105059276.FNA","EGU105059342.FNA","EGU105059366.FNA",
"EGU105059367.FNA","EGU105059381.FNA","EGU105059441.FNA","EGU105059453.FNA","EGU105059458.FNA","EGU105059479.FNA","EGU105059480.FNA","EGU105059490.FNA",
"EGU105059570.FNA","EGU105059573.FNA","EGU105059575.FNA","EGU105059587.FNA","EGU105059594.FNA","EGU105059612.FNA","EGU105059624.FNA","EGU105059636.FNA",
"EGU105059639.FNA","EGU105059671.FNA","EGU105059853.FNA","EGU105059900.FNA","EGU105059996.FNA","EGU105060095.FNA","EGU105060589.FNA","EGU105061025.FNA","EGU105061385.FNA",
"EGU105061427.FNA","HEY1007.FNA","HEY1013.FNA","HEY1017.FNA","HEY1020.FNA","HEY1025.FNA","HEY1035.FNA","HEY1050.FNA",
"HEY1052.FNA","HEY1064.FNA","HEY110.FNA","HEY1119.FNA","HEY1168.FNA","HEY1171.FNA","HEY1197.FNA","HEY1201.FNA",
"HEY120.FNA","HEY122.FNA","HEY125.FNA","HEY12.FNA","HEY136.FNA","HEY139.FNA","HEY1484.FNA","HEY148.FNA",
"HEY1494.FNA","HEY14.FNA","HEY150.FNA","HEY1615.FNA","HEY164.FNA","HEY168.FNA","HEY17.FNA","HEY1801.FNA",
"HEY180.FNA","HEY1815.FNA","HEY182.FNA","HEY1842.FNA","HEY1854.FNA","HEY1877.FNA","HEY1901.FNA","HEY191.FNA",
"HEY194.FNA","HEY197.FNA","HEY1986.FNA","HEY201.FNA","HEY204e.FNA","HEY204s.FNA","HEY2056.FNA","HEY207.FNA",
"HEY215.FNA","HEY2164.FNA","HEY218.FNA","HEY21.FNA","HEY2238.FNA","HEY225.FNA","HEY2291.FNA","HEY231.FNA",
"HEY2339.FNA","HEY2363.FNA","HEY2370.FNA","HEY2377.FNA","HEY237.FNA","HEY2388.FNA","HEY240.FNA","HEY2459.FNA",
"HEY245.FNA","HEY24.FNA","HEY250.FNA","HEY252e.FNA","HEY252p.FNA","HEY252s.FNA","HEY2550.FNA","HEY2561.FNA",
"HEY257.FNA","HEY267.FNA","HEY269.FNA","HEY277.FNA","HEY280.FNA","HEY281.FNA","HEY282.FNA","HEY290.FNA",
"HEY293.FNA","HEY296.FNA","HEY299.FNA","HEY305.FNA","HEY308.FNA","HEY310.FNA","HEY31.FNA","HEY323.FNA",
"HEY326.FNA","HEY32e.FNA","HEY32s.FNA","HEY332.FNA","HEY334.FNA","HEY340.FNA","HEY357.FNA","HEY360.FNA",
"HEY362.FNA","HEY363.FNA","HEY369.FNA","HEY378e.FNA","HEY378s.FNA","HEY38.FNA","HEY391.FNA","HEY392.FNA",
"HEY415.FNA","HEY417.FNA","HEY421.FNA","HEY429.FNA","HEY449.FNA","HEY464.FNA","HEY484.FNA","HEY490.FNA",
"HEY497.FNA","HEY4.FNA","HEY508.FNA","HEY514.FNA","HEY51.FNA","HEY52.FNA","HEY556.FNA","HEY563.FNA",
"HEY576.FNA","HEY581.FNA","HEY587.FNA","HEY604.FNA","HEY609.FNA","HEY61.FNA","HEY629.FNA","HEY630.FNA",
"HEY637.FNA","HEY673.FNA","HEY680.FNA","HEY703.FNA","HEY717.FNA","HEY727.FNA","HEY728.FNA","HEY732.FNA",
"HEY736.FNA","HEY740.FNA","HEY743.FNA","HEY757.FNA","HEY758.FNA","HEY762.FNA","HEY763.FNA","HEY785.FNA",
"HEY790.FNA","HEY793.FNA","HEY7.FNA","HEY807.FNA","HEY808.FNA","HEY822.FNA","HEY825.FNA","HEY82.FNA",
"HEY83.FNA","HEY84.FNA","HEY855.FNA","HEY856.FNA","HEY863.FNA","HEY872.FNA","HEY874.FNA","HEY883e.FNA",
"HEY883n.FNA","HEY886.FNA","HEY88.FNA","HEY897.FNA","HEY89.FNA","HEY938.FNA","HEY948.FNA","HEY94.FNA",
"HEY950.FNA","HEY958.FNA","HEY964.FNA","HEY977.FNA","HEY982.FNA","HEY985.FNA","HEY989.FNA",]

## Hybpiper with intronerate function
for i in range(0, len(sp)):
    gwf.target_from_template('Hybpiper_'+str(i), hybpiper(sp = sp[i],
                                                        p1 = "_clean-READ1.fastq",
                                                        p2 = "_clean-READ2.fastq",
                                                        un = "_clean-READ12-single.fastq",
                                                        path_out= "/home/kris/bp_ceroxyloideae/hybpiper/",
                                                        path_in = "/home/kris/palms_phylogeny/5.Ceroxyloids/",
                                                        done = "/home/kris/bp_ceroxyloideae/done_hybpiper/"+sp[i]))


## Generating stats
gwf.target_from_template('stats', stats(path_in = "/home/kris/bp_ceroxyloideae/hybpiper/"))



## Generating Heatmap
gwf.target_from_template('heatmap', heatmap(path_in = "/home/kris/bp_ceroxyloideae/stats/"))


## Generating Paralogs 2
gwf.target_from_template('paralogs2', paralogs2(path_in = "/home/kris/bp_ceroxyloideae/hybpiper/"))



#### Coverage
for i in range(0, len(sp)):
    gwf.target_from_template('Coverage2_'+sp[i], coverage2(sp2 = sp[i],
                                                        path_in = "/home/kris/bp_ceroxyloideae/hybpiper/",
                                                        all_bam = "_all.bam",
                                                        all_sorted_bam ="_all_sorted.bam",
                                                        all_sorted_bam_bai="_all_sorted.bam.bai",
                                                        bam =".bam",
                                                        cov=".cov",
                                                        fasta = ".fasta",
                                                        fasta_amb = ".fasta.amb",
                                                        fasta_ann = ".fasta.ann",
                                                        fasta_bwt = ".fasta.bwt",
                                                        fasta_pac = ".fasta.pac",
                                                        fasta_sa = ".fasta.sa",
                                                        trimmed_fasta = "_trimmed.fasta",
                                                        up_bam = "_up.bam",
                                                        path_out = "/home/kris/bp_ceroxyloideae/coverage/",
                                                        done = "/home/kris/bp_ceroxyloideae/coverage/done_coverage/"+sp[i],
                                                        dir_in ="/home/kris/palms_phylogeny/5.Ceroxyloids/", #Folder with clean reads + unpaired
                                                        dir_out ="/home/kris/bp_ceroxyloideae/coverage/", # folder with coverage
                                                        dir_wrk = "/home/kris/bp_ceroxyloideae/hybpiper/" ))



#### Retrieve sequences and sort into files with gene names 2
gwf.target_from_template('retrieve2', retrieve2(path_in="/home/kris/bp_ceroxyloideae/coverage/"))


#### MAFFT 2
for i in range(len(genes)):
    gwf.target_from_template('Mafft_'+str(i), mafft(genes = genes[i],
                                                        path_out= "/home/kris/bp_ceroxyloideae/mafft/",
                                                        path_in = "/home/kris/bp_ceroxyloideae/retrieve/",
                                                        done = "/home/kris/bp_ceroxyloideae/done_mafft/"+genes[i]))

#### Trimal according to a pre-defined gt values
for i in range(len(genes)):
    gwf.target_from_template('gt_trimming_'+genes[i], gt_trimming(genes = genes[i],
                                                        path_in = "/home/kris/bp_ceroxyloideae/mafft/",
                                                        done = "/home/kris/bp_ceroxyloideae/done_trimal/"+genes[i]))

### AMAS statistics
gwf.target_from_template('amas_raw', amas_raw(path_in="/home/kris/bp_ceroxyloideae/trimal/"))



#### Optrimal
gwf.target_from_template('optrim', optrim(path_in="/home/kris/bp_ceroxyloideae/trimal/"))


###Running CIAlign on the trimmed_fasta - Including Paralogs
for i in range(0, len(genes)):
    gwf.target_from_template('Cialign1'+str(i), cialign1(genes = genes[i],
                                            path_in = "/home/kris/bp_ceroxyloideae/best_alignments/optimal_final_results/",
                                            path_out = "/home/kris/bp_ceroxyloideae/cialign/"))

## Running TAPER after CIALIGN
for i in range(0, len(genes)):
    gwf.target_from_template('Taper_'+str(i), taper(genes = genes[i],
                                                    path_in = "/home/kris/bp_ceroxyloideae/cialign/",
                                                    path_out = "/home/kris/bp_ceroxyloideae/taper/"))


#Running IQTREE for files trimmed with trimal and CIAlign                                             
for i in range(0, len(genes)):
   gwf.target_from_template('Iqtree_'+str(i), iqtree(genes = genes[i],

                                                    path_in = "/home/kris/bp_ceroxyloideae/taper/"))


# Running ASTRAL f
gwf.target_from_template('astral_2', astral_2(genes = genes[i],
                                                    path_in = "/home/kris/bp_ceroxyloideae/astral/",
                                                    gene_tree_file="gene_trees.nex",
                                                    output="astral_tree_probabilities.tre"))
                                                        
# Running ASTRAL f
gwf.target_from_template('astral_gene_supporting', astral_gene_supporting(genes = genes[i],
                                                    path_in = "/home/kris/bp_ceroxyloideae/astral/",
                                                    gene_tree_file="genes_trees.nex",
                                                    output="astral_gene_probability.tre"))
                                                        


