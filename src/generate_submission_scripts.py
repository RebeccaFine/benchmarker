import os, sys
import ConfigParser

# input gene list name must be formatted as: TRAITNAME_LABEL_Genes.txt
# one label consistently used for each gene file

Config = ConfigParser.ConfigParser()
Config.read(sys.argv[1])

output_directory = Config.get("INPUT_FILES_AND_SETTINGS",'output_directory')
if not output_directory.endswith('/'):
    output_directory = output_directory + '/'

annotation_script_path = Config.get("SCRIPT_PATHS",'annotation_script_path')
ldsc_script_path = Config.get("SCRIPT_PATHS", 'ldsc_script_path')
list_of_gwas_paths = Config.get( "INPUT_FILES_AND_SETTINGS", 'list_of_gwas_paths_and_traits')
label = Config.get( "INPUT_FILES_AND_SETTINGS", 'label')
prio_gene_list  = Config.get( "INPUT_FILES_AND_SETTINGS", 'prio_gene_list')
prio_gene_list_has_header= Config.getboolean("INPUT_FILES_AND_SETTINGS","prio_gene_list_has_header")
gene_col_name = Config.get("INPUT_FILES_AND_SETTINGS", 'gene_col_name')
if gene_col_name == '':
    gene_col_name = 'NA'
window_size_in_kb = Config.get( "INPUT_FILES_AND_SETTINGS", 'window_size_in_kb')
bfile_path = Config.get("DATA_FILES", 'bfile_path')
print_snps_path = Config.get("DATA_FILES", 'print_snps_path')
frq_path = Config.get("DATA_FILES", 'frq_path')
weights_path = Config.get("DATA_FILES", 'weights_path')
baseline_path = Config.get("DATA_FILES", 'baseline_path')
gene_boundary_path = Config.get("DATA_FILES", 'gene_boundary_path')
gene_window_path = Config.get( "INPUT_FILES_AND_SETTINGS",'gene_window_path')
normalize_script_path = Config.get("SCRIPT_PATHS",'normalize_script_path')
concatenation_script_path = Config.get("SCRIPT_PATHS", 'concatenation_script_path')
template_file_annotate = Config.get('SCRIPT_PATHS', 'template_file_annotate')
template_file_partition = Config.get('SCRIPT_PATHS','template_file_partition')
main_script_directory = Config.get("SCRIPT_PATHS", 'main_script_directory')
if not main_script_directory.endswith('/'):
    main_script_directory = main_script_directory + '/'

template_file_annotate =  main_script_directory + template_file_annotate
template_file_partition =  main_script_directory + template_file_partition 
annotation_script_path =  main_script_directory + annotation_script_path
normalize_script_path  =  main_script_directory + normalize_script_path 
concatenation_script_path = main_script_directory + concatenation_script_path

gwas_paths = []
print list_of_gwas_paths
with open(list_of_gwas_paths) as x:
    x.readline()
    for line in x:
        gwas_path, trait = line.strip().split('\t')
        gwas_paths.append( (gwas_path, trait) )

joblist = '{label}_jobList.txt'.format( label = label)
with open( joblist, 'w') as x:
    x.write('gwas_path' + '\t' + 'trait' + '\t' + 'chrom' + '\n')
    for gwas_path, trait in gwas_paths:
        for chrom in range(1,23):
            x.write( gwas_path + '\t' + trait + '\t' + str(chrom)  + '\n')

joblist_noChrom = '{label}_jobList_noChrom.txt'.format( label = label)
with open( joblist_noChrom, 'w') as x:
    x.write('gwas_path' + '\t' + 'trait' '\n')
    for gwas_path, trait in gwas_paths:
        x.write( gwas_path + '\t' + trait + '\n')
    
qsub_script_1 = '{output_directory}/{label}_annotateAndComputeLDScores.sh'.format( label = label, output_directory = output_directory)
qsub_script_2 = '{output_directory}/{label}_partitionHeritability.sh'.format( label = label,output_directory = output_directory)
script_3 = '{output_directory}/{label}_concatenate.sh'.format(label = label,output_directory = output_directory)

print 'scripts:'
print qsub_script_1 
print qsub_script_2
print script_3

with open( qsub_script_1, 'w') as outfile:
    with open( template_file_annotate ) as template:
        for line in template:
            if line.startswith('#$ -t'):
                outline = '#$ -t 1-{job_num}'.format( job_num = len(gwas_paths) * 22 ) 
            elif line.startswith('output_directory'):
                outline = 'output_directory={output_directory}'.format( output_directory = output_directory)
            elif line.startswith( 'annotate_script_path=' ):
                outline = 'annotate_script_path={annotation_script_path}'.format(annotation_script_path = annotation_script_path)
            elif line.startswith('ldsc_script_path'):
                outline = 'ldsc_script_path={ldsc_script_path}'.format(ldsc_script_path = ldsc_script_path)
            
            elif line.startswith('joblist'):
                outline = 'joblist={joblist}'.format(joblist = joblist)
           
            elif line.strip() == 'prio_gene_list=':
                outline = 'prio_gene_list={prio_gene_list}'.format( prio_gene_list = prio_gene_list)
            elif line.startswith('baseline_path'):
                outline = 'baseline_path=%s.${curr_chromosome}.' %(baseline_path.rstrip('.'))
            elif line.startswith('gene_boundary_path'):
                outline = 'gene_boundary_path={gene_boundary_path}'.format( gene_boundary_path = gene_boundary_path )
            elif line.startswith('bfile_path'):
                outline = 'bfile_path=%s.${curr_chromosome}' %(bfile_path.rstrip('.'))
              
            elif line.startswith('label'):
                outline = 'label={label}'.format( label = label )
            elif line.startswith('window_size'):
                outline = 'window_size_in_kb={window_size_in_kb}'.format( window_size_in_kb = window_size_in_kb)
            elif line.startswith('print_snps_path'):
                outline = 'print_snps_path={print_snps_path}'.format( print_snps_path = print_snps_path)
            elif line.strip() == 'prio_gene_list_has_header=':
                outline = 'prio_gene_list_has_header={prio_gene_list_has_header}'.format( prio_gene_list_has_header = prio_gene_list_has_header)
            elif line.startswith('gene_col_name'):
                outline = 'gene_col_name={gene_col_name}'.format( gene_col_name  = gene_col_name )
            else:
                outline = line.strip()

            outfile.write(outline + '\n')

with open(qsub_script_2, 'w') as outfile:
       with open( template_file_partition ) as template:
        for line in template:
            if line.startswith('#$ -t'):
                outline = '#$ -t 1-{job_num}'.format( job_num = len(gwas_paths)  )
            elif line.startswith('ldsc_script_path'):
                outline = 'ldsc_script_path={ldsc_script_path}'.format(ldsc_script_path = ldsc_script_path)
    
            elif line.startswith('joblist'):
                outline = 'joblist={joblist_noChrom}'.format(joblist_noChrom = joblist_noChrom)
            elif line.startswith('output_directory'):
                outline = 'output_directory={output_directory}'.format( output_directory = output_directory)
    
            elif line.startswith('baseline_path'):
                outline = 'baseline_path={baseline_path}'.format( baseline_path = baseline_path)
            elif line.startswith('weights_path'):
                outline = 'weights_path={weights_path}'.format( weights_path = weights_path) 
            elif line.startswith('label'):
                outline = 'label={label}'.format( label = label )
            elif line.startswith('window_size'):
                outline = 'window_size_in_kb={window_size_in_kb}'.format( window_size_in_kb = window_size_in_kb)
            elif line.startswith('frq_path'):
                outline = 'frq_path={frq_path}'.format( frq_path = frq_path )
            elif line.startswith('gene_window_path'):
                outline = 'gene_window_path={gene_window_path}'.format( gene_window_path = gene_window_path )
            elif line.startswith('normalize_script_path='):
                outline = 'normalize_script_path={normalize_script_path} \\'.format (normalize_script_path = normalize_script_path )
            else:
                outline = line.strip()

            outfile.write(outline + '\n') 

with open( script_3, 'w') as outfile:
    outline = 'python {concatenation_script_path} \\'.format( concatenation_script_path = concatenation_script_path)
    outfile.write(outline + '\n')   
    outline = '--label {label} \\'.format( label = label )
    outfile.write(outline + '\n')
    outline = '--joblist {joblist} \\'.format( joblist = joblist_noChrom)
    outfile.write(outline + '\n')
    outline = '--output_directory {output_directory}'.format( output_directory = output_directory)
    outfile.write(outline + '\n')

cmd = 'chmod +x {script_3}'.format( script_3 = script_3)
os.system(cmd)
