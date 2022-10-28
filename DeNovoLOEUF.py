import hail as hl
hl.init()

if not hl.hadoop_exists('cohort.mt'): 
    mt = hl.read_matrix_table('cohortt.mt', _n_partitions=32)
    mt.write('cohort_32partitions.mt')
mt = hl.read_matrix_table('cohort_32partitions.mt')

#reduce dataset down
mt = mt.filter_rows(mt.gnomad_exomes.AF <= 0.001)

#load in gnomAD frequencies
#AF QC
mt = hl.split_multi_hts(mt)
mt = hl.variant_qc(mt)

#load in pedigree files
pedigree = hl.Pedigree.read('data/trio.ped') #pedigree data of callset

#generate de novo callset for later
de_novo_callset = hl.de_novo(mt, pedigree, mt.variant_qc.AF[1])
# de_novo_callset.describe()

#make abridged trio dataset
trio_abr = hl.trio_matrix(mt, pedigree, complete_trios=True)
# trio_abr.describe()

#import disease genes needed for later
if not hl.hadoop_exists('genes.ht'): 
    genes = hl.import_table('disease_genes.tsv', impute=True)
    genes = genes.key_by('geneIds')
    genes.write('genes.ht', overwrite=True)
genes = hl.read_table('genes.ht')
genes = genes.key_by('geneIds')


#import families and their clinical data
if not hl.hadoop_exists('CMG.ht'):
    CMG = hl.import_table('clinical.tsv', impute=True)
    CMG.key_by('s')
    CMG.write('CMG.ht')
CMG = hl.read_table('CMG.ht')

#now add in family data
mt = mt.annotate_cols(cmg_data = CMG.key_by('s')[mt.s])

#now add in GenCC data
rows=mt.rows()
gnomad = rows.select('gnomad_exomes')
rows = rows.select('geneIds')
rows = rows.explode('geneIds')
rows = rows.annotate(GenCC = genes[rows.geneIds])

mt = mt.annotate_rows(
    GenCC = rows.index(mt.row_key, all_matches=True).GenCC
)

#make trio unabridged dataset with key info needed for tool
trio = hl.trio_matrix(mt, pedigree, complete_trios=True)

#######################
#CODING VARIANT FILTER#
#######################

# filter callset first to max frequency needed
trio2=trio.filter_rows((trio.gnomad_exomes.AF < 0.001))

#filter for dom disease gene + AF < 0.001 + pLoF with LOEUF <0.2

coding_dom=trio2.filter_rows(

    (
        ((trio2.GenCC[0].inheritance_class == 'Dom') & (trio2.gnomad_exomes.AF < 0.001)) &
        (
            (
                (trio2.GenCC[0].LOEUF <0.2) & (trio2.mainTranscript.category == 'lof')) |
            (
                (
                    (trio2.mainTranscript.major_consequence == 'splice_region_variant') |
                    (trio2.mainTranscript.major_consequence == 'splice_donor_variant') |
                    (trio2.mainTranscript.major_consequence == 'splice_acceptor_variant')) &
                (trio2.splice_ai.delta_score > 0.2)) )))


#select probands only
filter_condition = (coding_dom.proband_entry.GT.is_non_ref())

coding_dom_proband = coding_dom.filter_entries(filter_condition)

#remove lines where all probands are NA
coding_dom_proband = coding_dom_proband.filter_rows(
    hl.agg.all(
        hl.is_missing(coding_dom_proband.proband_entry.GT)
    ),keep=False
)

#Genotypes for all individuals with parents
coding_dom_proband=coding_dom_proband.explode_rows(coding_dom_proband.genotypes)
