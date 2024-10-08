import gffutils

db = gffutils.create_db("/ref/annotation/hg38.ncbiRefSeq.gtf",
                        dbfn="/ref/qc/rseqc/annotation.db",
                        force=True,
                        keep_order=True,
                        merge_strategy='merge',
                        sort_attribute_values=True,
                        disable_infer_genes=True,
                        disable_infer_transcripts=True)

with open("/work/qc/rseqc/annotation.bed", 'w') as outfileobj:
    for tx in db.features_of_type('transcript', order_by='start'):
        bed = [s.strip() for s in db.bed12(tx).split('\t')]
        bed[3] = tx.id
        outfileobj.write('{}\n'.format('\t'.join(bed)))
