import os
import glob
import shutil
import lsst.sims.maf.batches as batches
import lsst.sims.maf.db as db
import lsst.sims.maf.metricBundles as mb
import argparse


if __name__ == "__main__":
    """
    Run the glance batch on all .db files in a directory.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=str, default=None)
    args = parser.parse_args()

    if args.db is None:
        if os.path.isfile('trackingDb_sqlite.db'):
            os.remove('trackingDb_sqlite.db')
        db_files = glob.glob('*.db')
    else:
        db_files = [args.db]
    run_names = [os.path.basename(name).replace('.db', '') for name in db_files]

    for filename, name in zip(db_files, run_names):
        if os.path.isdir(name):
            shutil.rmtree(name)
        opsdb = db.OpsimDatabaseV4(filename)
        colmap = batches.ColMapDict('OpsimV4')

        bdict = {}
        bdict.update(batches.glanceBatch(colmap, name))
        bdict.update(batches.fOBatch(colmap, name))
        resultsDb = db.ResultsDb(outDir=name)
        group = mb.MetricBundleGroup(bdict, opsdb, outDir=name, resultsDb=resultsDb, saveEarly=False)
        group.runAll(clearMemory=True, plotNow=True)
        resultsDb.close()
        opsdb.close()
        db.addRunToDatabase(name, 'trackingDb_sqlite.db', None, name, '', '', name+'.db')
