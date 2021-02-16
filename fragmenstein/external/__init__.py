import importlib
import json
import os

import sys
from typing import List, Dict


class ExternalToolImporter():

  #Default config file will be used if --external_tools_config_file flag is not provided
  if "--external_tools_config_file" in sys.argv: #TODO; change it for an environmental variable
    configFname = sys.argv[sys.argv.index("--external_tools_config_file") + 1]
    configFname = os.path.abspath(os.path.expanduser(configFname))
  else:
    configFname = os.path.join(os.path.split(__file__)[0], "external_config.json") #Modify this attribute if using another file

  config_json_fname = configFname

  with open(config_json_fname) as f:
    config_data = json.load( f )

  @classmethod
  def _get_data_dict(cls, toolName: str) -> Dict[str, Dict]:
    '''
    :param toolName: str. The external tool name
    :return: the config dict associated to toolName
    '''
    try:
      return ExternalToolImporter.config_data[toolName]
    except KeyError:
      raise  ValueError("Error, desired tool is not included in config file (%s)" % ExternalToolImporter.config_json_fname)

  @classmethod
  def get_rootdir(cls, toolName):
    '''
    :param toolName: str. The external tool name
    :return: the ROOT_DIR for the toolName
    '''
    return  os.path.expanduser( cls._get_data_dict(toolName)["ROOT_DIR"] )

  @classmethod
  def import_tool(cls, toolName: str, moduleNames: List[str]) -> List['class_module']:
    '''
    :param toolName: str. The external tool name
    :param moduleNames: the modules that want to be loaded from toolName. They can be submodules. E.g. externalName1.tools.foo
    :return: a list of modules loaded

    '''
    data = cls._get_data_dict(toolName)
    data["ROOT_DIR"] = os.path.expanduser(data["ROOT_DIR"] )

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
  print("this is a test")
  ExternalToolImporter.import_tool("DeLinker", ["DeLinker_test", "data.prepare_data"])
  [pyrosetta]= ExternalToolImporter.import_tool("pyrosetta", ["pyrosetta"])
  print( pyrosetta , type(pyrosetta))