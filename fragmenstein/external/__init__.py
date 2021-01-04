import importlib
import json
import os

import sys
from typing import List, Dict


class ExternalToolImporter():

  config_json_fname = os.path.join(os.path.split(__file__)[0], "external_config.json") #Modify this attribute if using another file

  with open(config_json_fname) as f:
    config_data = json.load( f )


  @classmethod
  def _get_data_dict(cls, toolName: str) -> Dict[str, Dict]:
    try:
      return ExternalToolImporter.config_data[toolName]
    except KeyError:
      raise  ValueError("Error, desired tool is not included in config file (%s)" % ExternalToolImporter.config_json_fname)

  @classmethod
  def get_rootdir(cls, toolName):
    return  cls._get_data_dict(toolName)["ROOT_DIR"]

  @classmethod
  def import_tool(cls, toolName: str, moduleNames: List[str]) -> List['class_module']:
    data = cls._get_data_dict(toolName)
    for dir in data["PYTHONPATH"]:
      dir = dir % data
      if dir not in sys.path:
        sys.path.append(dir)
    modules=[]
    for moduleName in moduleNames:
      modules.append(importlib.import_module(moduleName))

    for code in data["POST_IMPORT"]:
      exec(code)
    return modules

if __name__=="__main__":
  ExternalToolImporter.import_tool("DeLinker", ["DeLinker_test"])
  [pyrosetta]= ExternalToolImporter.import_tool("pyrosetta", ["pyrosetta"])
  print( pyrosetta , type(pyrosetta))