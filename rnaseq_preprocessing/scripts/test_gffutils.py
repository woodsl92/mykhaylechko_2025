import tempfile
import gffutils
import os
from argparse import ArgumentParser

def get_gtf_db(gtf, in_memory=False):
    """
    create a gffutils DB
    """
    db_file = gtf + ".db"
    if os.path.exists(db_file):
        return gffutils.FeatureDB(db_file)
    db_file = ":memory:" if in_memory else db_file
    if in_memory or not os.path.exists(db_file):
        infer_extent = guess_infer_extent(gtf)
        db = gffutils.create_db(gtf, dbfn=db_file,
                                infer_gene_extent=infer_extent)
    if in_memory:
        return db
    else:
        return gffutils.FeatureDB(db_file)

def guess_infer_extent(gtf_file):
    """
    guess if we need to use the gene extent option when making a gffutils
    database by making a tiny database of 1000 lines from the original
    GTF and looking for all of the features
    """
    _, ext = os.path.splitext(gtf_file)
    tmp_out = tempfile.NamedTemporaryFile(suffix=".gtf", delete=False).name
    with open(tmp_out, "w") as out_handle:
        count = 0
        in_handle = open(gtf_file) if ext != ".gz" else gzip.open(gtf_file)
        for line in in_handle:
            if count > 1000:
                break
            out_handle.write(line)
            count += 1
        in_handle.close()
    db = gffutils.create_db(tmp_out, dbfn=":memory:", infer_gene_extent=False)
    os.remove(tmp_out)
    features = [x for x in db.featuretypes()]
    if "gene" in features and "transcript" in features:
        return False
    else:
        return True

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("gtf_file")
    args = parser.parse_args()
    db = get_gtf_db(args.gtf_file)
    for x in db.featuretypes():
        print(x)