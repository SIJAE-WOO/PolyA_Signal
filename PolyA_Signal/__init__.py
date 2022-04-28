import os
import tabix

DeepPASTA_path = '/home/sijaewoo/processed_data/polyA_site_prediction/genome_wide_polyA_site_prediction_human.hg38.sorted.bed.gz'
APADB_path = '/home/sijaewoo/processed_data/polyA_site_prediction/apadb_v2_final.hg38.sorted.bed.gz'
    

def PolyA_Signal(loc, data='DeepPASTA'):

    chrN,tok = loc.split(':')
    posSta,posEnd = list(map(int,tok.split('-')))

    if data == 'DeepPASTA':
        polya_path = DeepPASTA_path
        tb_polya = tabix.open(polya_path)
        result_records = tb_polya.querys('chr%s:%s-%s' % (chrN,posSta-100,posEnd+100))
    elif data == 'APADB':
        polya_path = APADB_path
        tb_polya = tabix.open(polya_path)
        result_records = tb_polya.querys('chr%s:%s-%s' % (chrN,posSta,posEnd+1))

    resultL = []
    polyAL = []

    for result in result_records:
        resultL.append((result[:]))
        
    if data == 'DeepPASTA':
        for i in resultL:
            polyA_loc = ('%s:%s-%s' % (i[0], i[1], i[2]))
            polyAL.append((polyA_loc, i[3], i[4]))
    
    elif data == 'APADB':
        for i in resultL:
            polyA_loc = ('%s:%s-%s' % (i[0], i[1], i[2]))
            gene = i[3].split('.')[0]
            region = i[3].split(':')[1]
            polyA_no = i[3].split('.')[1].split(':')[0]
            polyAL.append((polyA_loc, i[5], gene, polyA_no, region))
        
    return polyAL