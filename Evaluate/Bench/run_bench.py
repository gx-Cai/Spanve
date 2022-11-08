import os
import tqdm

def filter_datasets(name):
    if not filter_datasets_flag: 
        if filter_code in name:
            return True
        else: return False


filter_datasets_flag = False
check_runned = True
filter_code = ''

saving_root = '../recode'
datasets_id = '10XGenomics'
base_dir = f'../Data/{datasets_id}'

dataset_dirs = [os.path.join(base_dir,i) for i in os.listdir(base_dir) if filter_datasets(i)]

softs = [] 
## all Available softs:
# -- Run on python -- 
# 'p.spanve-k','p.spanve-d' # spanve with preprocess
# 'spanve-k','spanve-d','sepal','somde','moran','geary','spatialde'
# -- Run on R --
# 'r-sparkx','r-gitto-rank','r-gitto-sirank','r-meringue'

bar = tqdm.tqdm(total=len(dataset_dirs)*len(softs))

for dataset_dir in dataset_dirs:
    if not dataset_dir.endswith('h5ad'): 
        bar.update(1)
        continue
    for soft in softs:
        bar.set_description(f'{dataset_dir}_{soft}')
        dataset_name = os.path.basename(dataset_dir)
        saving_path = os.path.join(saving_root,datasets_id,dataset_name,soft)
        if not os.path.exists(saving_path):
            os.makedirs(saving_path)
        elif check_runned:
            bar.update(1)
            continue

        memory_info = os.path.join(saving_path,'memory_info.txt')
        command = f'\\time --output {memory_info} python ./frame.py --dataset {dataset_dir} --soft {soft} --saving_path {saving_path}'
        print(command)     
        try:
            os.system(command)
        except:
            raise Exception(f'error in {dataset_dir} {soft}')
            os.exit(1)
        bar.update(1)
bar.close()