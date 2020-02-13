import glob

# To generate sym links:
# find . -name *10yrs.db | xargs -I'{}' ln -s '{}' maf_all/.

if __name__ == "__main__":
    db_files = glob.glob('*10yrs.db')
    outfile = open("run_all.sh", 'w')

    for filename in db_files:
        outfile.write('python ../glance_dir.py --db %s \n' % filename)
        outfile.write('python ../scimaf_dir.py --db %s \n' % filename)

    outfile.close()
