# https://stackoverflow.com/questions/19425736/how-to-redirect-stdout-and-stderr-to-logger-in-python

class LoggerWriter:
    """
    Write stderr to logger.
    See capture_rdkit_log method of Victor
    where this happens:
    sys.stderr = LoggerWriter(cls.journal.warning)
    """
    def __init__(self, writer):
        self._writer = writer
        self._msg = ''

    def write(self, message):
        self._msg += str(message)
        # this causes an infinite recursion and
        # I am not sure why it was added
        # while '\n' in self._msg:
        #     pos = self._msg.find('\n')
        #     self._writer(self._msg[:pos])
        #     self._msg = self._msg[pos + 1:]

    def flush(self):
        if self._msg != '':
            self._writer(self._msg)
            self._msg = ''
