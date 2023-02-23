from typing import *
from ._victor_base import _VictorBase
import warnings, time

class _VictorSafety(_VictorBase):

    # cvar: error_to_catch = Exception

    def _safely_do(self,
                   execute: Optional[Callable] = None,
                   resolve: Optional[Callable] = None,
                   reject: Optional[Callable] = None):
        """
        A safety net around the analysis.
        Ought to be a decorator and ought to not use the same names as a JS Promise.
        The class attribute ``error_to_catch`` is by default Exception

        :param execute: what to run (main)
        :param resolve: what to run at the end (regardless of failure)
        :param reject: what to run if ``exceute`` fails
        :return:
        """
        # warnings
        with warnings.catch_warnings(record=True) as self._warned:
            try:
                if execute is not None:
                    execute()
            except KeyboardInterrupt as err:
                raise err
            except self.error_to_catch as err:
                self.error_msg = f'{err.__class__.__name__}: {err}'
                if reject is not None:
                    reject(err)
            finally:
                if resolve is not None:
                    resolve()

    def _resolve(self) -> None:
        """
        This gets called at the end of ``_safely_do``, regardless of the whether it worked or not.
        So the name is a bit misleading.

        :return:
        """
        self.tock = time.time()
        self.journal.info(f'{self.long_name} - Time taken: {self.tock - self.tick}')

    def _reject(self, err) -> None:
        """
        This gets called by ``_safely_do`` on error.

        :param err: the error raised.
        :return:
        """
        self.journal.error(f'{self.long_name} — {err.__class__.__name__}: {err}')