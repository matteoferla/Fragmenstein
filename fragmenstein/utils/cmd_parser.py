import os
import argparse

from fragmenstein.utils.config_manager import ConfigManager


class ArgumentParser(argparse.ArgumentParser):
    '''
    A subclass of argparse.ArgumentParser that automatically adds arguments defined in ConfigManager and
    saves the provided values for those default values in os.environ, so that ConfigManager get can retrieve them
    '''

    @classmethod
    def from_nameInConfig_to_argname(cls, name ):
        return "--"+name.lower()

    @classmethod
    def from_argname_to_nameInConfig(cls, name ):
        return name.strip("--").upper()

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)
        for (name, defaultVal), help in ConfigManager.yieldParseableParams():
            if defaultVal is None:
                varType = str
            else:
                varType = type(defaultVal)
            self.add_argument(self.from_nameInConfig_to_argname(name), type= varType , default= defaultVal, help=help)

    def parse_args(self, *args, **kwargs):

        parsed_args = super().parse_args(*args, **kwargs)
        for arg,val in vars(parsed_args).items():
            if self.from_argname_to_nameInConfig(arg) in ConfigManager.parseable_params:
                os.environ[self.from_argname_to_nameInConfig(arg)] = str(val)
        return parsed_args
