
id_convert = {}

with open('/Volumes/WorkingDrive_BWP/_bearNecropsyRNAseq_Sept2021/data/isoforms/z_old_preFeb2022/new_merge.combined.renamed.mapping.txt') as a:
    for line in a.readlines():
        line = line.rstrip().split(',')
        id_convert[line[0]] = line[1]

with open('/Users/blairperry/Downloads/new_merge.pbIDs.fa','w') as out:
    with open('/Users/blairperry/Downloads/new_merge.fa') as a:
        for line in a.readlines():
            if line[0] == '>':
                old_id = line.split(' ')[0][1:]
                new_id = '>' + id_convert[old_id]
                print >> out, new_id
            else:
                print >> out, line.rstrip()
