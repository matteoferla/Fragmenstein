from ._victor_safety import _VictorSafety
from ._loggerwriter import LoggerWriter
import logging, sys, os, re, requests, unicodedata
from rdkit import Chem
from rdkit_to_params import Params


class _VictorJournal(_VictorSafety):

    def _log_warnings(self):
        if len(self._warned):
            for w in self._warned:
                self.journal.warning(f'{self.long_name} - {w.message} ({w.category})')
            self._warned.clear()

    @classmethod
    def enable_stdout(cls,
                      level=logging.INFO,
                      captured: bool = True) -> None:
        """
        The ``cls.journal`` is output to the terminal.
        Running it twice can be used to change level.

        :param level: logging level
        :param captured: capture rdkit and pyrosetta?
        :return: None
        """
        cls.journal.handlers = [h for h in cls.journal.handlers if h.name != 'stdout']
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(level)
        handler.set_name('stdout')
        handler.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s'))
        cls.journal.addHandler(handler)
        # logging.getLogger('py.warnings').addHandler(handler)
        if captured:
            cls.capture_logs()

    @classmethod
    def enable_logfile(cls, filename='reanimation.log',
                       level=logging.INFO,
                       captured: bool = True) -> None:
        """
        The journal is output to a file.
        Running it twice can be used to change level.

        :param filename: file to write.
        :param level: logging level
        :param captured: capture rdkit and pyrosetta?
        :return: None
        """
        cls.journal.handlers = [h for h in cls.journal.handlers if h.name != 'logfile']
        handler = logging.FileHandler(filename)
        handler.setLevel(level)
        handler.set_name('logfile')
        handler.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s'))
        cls.journal.addHandler(handler)
        # logging.getLogger('py.warnings').addHandler(handler)
        if captured:
            cls.capture_logs()

    _rdkit_captured = False

    @classmethod
    def capture_rdkit_log(cls):
        """
        RDKit spits a few warning and errors.
        This makes them inline with the logger.
        """
        if cls._rdkit_captured:
            return
        Chem.WrapLogs()
        sys.stderr = LoggerWriter(cls.journal.warning)
        cls._rdkit_captured = True

    _rosetta_captured = False

    @classmethod
    def capture_rosetta_log(cls):
        """
        Rosetta normally prints to stout. This captures the messages into ``journal``.
        It technically simply passes the handlers of journal to that of the Rosetta logger.
        For alternatives, https://github.com/matteoferla/pyrosetta_scripts/tree/main/init_helper
        """
        if cls._rosetta_captured:
            return
        import pyrosetta
        pyrosetta.logging_support.set_logging_sink()
        logger = logging.getLogger("rosetta")
        logger.setLevel(logging.DEBUG)
        logger.handlers = cls.journal.handlers
        cls._rosetta_captured = True

    @classmethod
    def capture_logs(cls):
        cls.capture_rdkit_log()
        cls.capture_rosetta_log()

    @classmethod
    def slack_me(cls, msg: str) -> bool:
        """
        Send message to a slack webhook

        :param msg: Can be dirty and unicode-y.
        :return: did it work?
        :rtype: bool
        """
        webhook = os.environ['SLACK_WEBHOOK']
        # sanitise.
        msg = unicodedata.normalize('NFKD', msg).encode('ascii', 'ignore').decode('ascii')
        msg = re.sub('[^\w\s\-.,;?!@#()\[\]]', '', msg)
        r = requests.post(url=webhook,
                          headers={'Content-type': 'application/json'},
                          data=f"{{'text': '{msg}'}}")
        if r.status_code == 200 and r.content == b'ok':
            return True
        else:
            return False


# ----------------------- Make params use the same log.

Params.log = _VictorJournal.journal
