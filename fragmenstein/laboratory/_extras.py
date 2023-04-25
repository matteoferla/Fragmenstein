class LabExtras:

    error_classifications ={'SUCCESS': 'success',
                            'UNCLASSIFIED': 'unclassified',
                            'ConnectionError': 'distance',  # legacy
                            'DistanceError': 'distance',
                            'RectificationError': 'incorrect rectification',
                            'FullOverlapError': 'distance',
                            'TimeoutError': 'timeout',
                            'No sub-structure match': 'incorrect parameterisation #1',
                            'Violation occurred': 'incorrect rectification #1',
                            'KekulizeException': 'incorrect rectification #1',
                            'AtomValenceException': 'incorrect rectification #1',
                            'Bad Conformer Id': 'incorrect rectification #1',
                            ('rosetta', 'not found'): 'incorrect parameterisation #2',
                            ('rosetta', 'nan'): 'incorrect parameterisation #3',
                            'but not the minimised': 'incorrect parameterisation #4',
                            'UtilityExitException': 'incorrect parameterisation #5',
                            'TypeError': 'embedding error',
                            'RecursionError': 'embedding error',
                            'utterly impossible': 'embedding error'
                      }
    @classmethod
    def error_classify(cls, text):
        if not text:
            return cls.error_classifications['SUCCESS']
        text = str(text)
        for key, value in cls.error_classifications.items():
            if isinstance(key, tuple):
                if all(x in text for x in key):
                    return value
            elif key in text:
                return value
            else:
                return cls.error_classifications['UNCLASSIFIED']

    def convert_to_sdf(self, df):
        pass
