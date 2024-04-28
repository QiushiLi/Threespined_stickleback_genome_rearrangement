import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'

gff = pd.read_csv("Gac-HiC_revised_genome_assembly.fa.tsv", sep="\t", comment='&', header=None)


gff['age_len'] = gff[6] - gff[5]


age_len_sum = gff.groupby(14)['age_len'].sum()

gff['age_mod'] = gff[15] * gff['age_len']
age_mod_sum = gff.groupby(14)['age_mod'].sum()

age_weighted = age_mod_sum / age_len_sum

start_min = gff.groupby(14)[5].min()
stop_max = gff.groupby(14)[6].max()

# Get unique scaffolds
scafs = gff.groupby(14)[4].unique()

# Process the split operation
gff['te_split'] = gff[9].apply(lambda x: x.split("#")[1] if len(x.split("#")) > 1 else None)

# Get min te_type per group
te_type = gff.groupby(14)['te_split'].min()

te = pd.DataFrame({
    'scafs': scafs,
    'start_min': start_min,
    'stop_max': stop_max,
    'te_type': te_type,
    'age_weighted': age_weighted
})

# Add the names of the scaffolds as a column
te['names_scafs'] = te.index
te.to_csv("temp_te.txt", sep="\t", index=False)
te = pd.read_csv("temp_te.txt", sep="\t", header=None)
stick_tube = pd.read_csv("stick_tube_mappings.txt", sep="\t", header=True)
not_moved = pd.read_csv("all_not_moved.txt", sep="\t", header=True)


not_moved.sort_values(by=[0, 1], inplace=True)  # Assuming 'stick_chrom' is the first column and 'stick_start' is the second
not_moved = not_moved[not_moved['tube_chrom'] != 'chrUN']

lsg = pd.read_csv("LSG_blast_completeness.txt", sep="\t", header=None)
res = pd.read_csv("rearrangement_blast_completeness.txt", sep="\t", header=None)

lsg = lsg[(lsg[3] > 0.9) | lsg[3].isna()]

res_not_te = pd.read_csv("rearranged_genes_notTE_stickle_denovo_TEmap0.1_no_chrUn.txt", sep="\t", header=True)
res_not_te = res_not_te[~res_not_te['stick_gene'].isin(res[0])]

missing_not_te = pd.read_csv("LSG_stickle_notTE_stickle_denovo_TEmap0.1_no_chrUn.txt", sep="\t", header=True)
dups_not_te = pd.read_csv("dups_notTE_stickle_denovo_TEmap0.1_no_chrUn_checkbass.txt", sep="\t", header=True)

temp1 = pd.read_csv("LSG_blast_rate.txt", sep="\t", header=None)
temp2 = pd.read_csv("ENSGACT_to_ENSGACP.txt", sep="\t", header=None)

temp2b = temp2[temp2[1].isin(temp1[0])]

missing_not_te['temp'] = missing_not_te['stick_gene'].apply(lambda x: x.split("Name=")[1] if "Name=" in x else None)
missing_not_te = missing_not_te[~(missing_not_te['temp'].isin(temp1[0]) | missing_not_te['temp'].isin(temp2b[0]))]
missing_not_te = missing_not_te[~missing_not_te['stick_gene'].isin(lsg[0])]

the_levels = sorted(missing_not_te['te_type'].unique())
the_chroms = not_moved['stick_chrom'].unique()


res1 = pd.DataFrame(index=range(10000), columns=range(7))

count2 = 1
stick_all = None

import numpy as np

stick_window = 2000000 

# the_chroms, te, not_moved, res_not_te, dups_not_te, missing_not_te

res1 = pd.DataFrame(columns=['chrom', 'pos', 'res_not_te_count', 'dups_not_te_count', 'missing_not_te_count'])
stick_all = pd.DataFrame()

for chrom in the_chroms:
    sub_stick = not_moved[not_moved['stick_chrom'] == chrom]
    sub_te = te[te['scafs'] == chrom]
    
    count1 = 0
    is_done = False
    
    while not is_done:
        # Select the relevant TE data within the current window
        which_te_stick = sub_te[(sub_te['start_min'] >= count1) & (sub_te['start_min'] <= (count1 + stick_window))]
        
        stick1 = which_te_stick['te_type'].value_counts().reindex(the_levels, fill_value=0)
        
        stick_all = stick_all.append(pd.DataFrame([[chrom, count1] + stick1.tolist()], columns=['chrom', 'pos'] + the_levels.tolist()))
        
        current_row = {
            'chrom': chrom,
            'pos': count1,
            'res_not_te_count': ((res_not_te['stick_chrom'] == chrom) & (res_not_te['stick_start'] >= count1) & (res_not_te['stick_start'] < (count1 + stick_window))).sum(),
            'dups_not_te_count': ((dups_not_te['stick_chrom'] == chrom) & (dups_not_te['stick_start'] >= count1) & (dups_not_te['stick_start'] < (count1 + stick_window))).sum(),
            'missing_not_te_count': ((missing_not_te['stick_chrom'] == chrom) & (missing_not_te['stick_start'] >= count1) & (missing_not_te['stick_start'] < (count1 + stick_window))).sum()
        }
        res1 = res1.append(current_row, ignore_index=True)
        
        if (count1 + stick_window) < sub_stick['stick_end'].max():
            count1 += stick_window
        else:
            is_done = True

stick_name = f"stick_simple_all_te_window_{stick_window // 1000}k.txt"
res1.dropna(subset=['chrom']).to_csv(stick_name, index=False, sep='\t', header=True, quotechar='"')

results_name = f"results_simple_all_te_window_{stick_window // 1000}k.txt"
res1.to_csv(results_name, index=False, sep='\t', header=True, quotechar='"')
