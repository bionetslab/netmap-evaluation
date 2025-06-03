
import yaml

def write_config(c, file):
    with open(file, "w") as handle:
        yaml.safe_dump(c, handle)


def split_index(aa):    
    aa.var['source']   = [l[0] for l in aa.var.index.str.split('_', expand=True)]
    aa.var['target']   = [l[1] for l in aa.var.index.str.split('_', expand=True)]
    return aa