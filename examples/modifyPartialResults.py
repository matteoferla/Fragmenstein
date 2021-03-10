import pickle

from fragmenstein.utils.io_utils import apply_func_to_files

#datadir = "/home/sanchezg/oxford/myProjects/diamondCovid/data/nsp13/enumeration/Site_1_brics/merges"
datadir = "/home/sanchezg/oxford/myProjects/diamondCovid/data/nsp13/enumeration/Site_7/brics/merges"

def fix_pickle(fname):
  with open(fname, "rb") as f:
    data = pickle.load(f)
  print(data)
  if isinstance(data, list):
    for comp in data:
      print(comp.ref_molIds)
      new_refMols =  sorted(set((x.primitiveId for x in comp.parents)))
      comp.ref_molIds = new_refMols
      print(comp.ref_molIds)
    with open(fname, "wb") as f:
      pickle.dump(data, f)

#data.ref_mols = [ "_".join(elem.split("_")[:2]) for elem in data.ref_mols]
#[x.primitiveId for x in data[0].parents]

apply_func_to_files(datadir, ".*\.final\.pickle", fix_pickle, use_parallel_dask= False)
