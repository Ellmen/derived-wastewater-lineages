import shutil

data_dir = '/Users/isaac/Code/covid/variant_analysis/data/fall_2021/ww_benchmark/samples/'

for i in range(100):
    s_dir = 'sample' + str(i+1)
    s_path = data_dir + s_dir
    # grom.coverage.csv grom.mapped.csv
    shutil.copyfile('{}/grom.coverage.csv'.format(s_path), './{}.coverage.csv'.format(s_dir))
    shutil.copyfile('{}/grom.mapped.csv'.format(s_path), './{}.mapped.csv'.format(s_dir))
