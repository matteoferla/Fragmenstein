'''

https://www.saltycrane.com/blog/2010/04/using-python-timeout-decorator-uploading-s3/

'''
import sys
import threading

from concurrent.futures import thread


def quit_function(fn_name):
    print('{0} took too long'.format(fn_name), file=sys.stderr)
    sys.stderr.flush() # Python 3 stderr is likely buffered.
    thread.interrupt_main() # raises KeyboardInterrupt


class TimeoutError(Exception):
    def __init__(self, value = "Timed Out"):
        self.value = value
    def __str__(self):
        return repr(self.value)

def timeout(seconds_before_timeout):
    def outer(fn):
        def inner(*args, **kwargs):
            timer = threading.Timer(seconds_before_timeout, quit_function, args=[fn.__name__])
            timer.start()
            try:
                result = fn(*args, **kwargs)
            finally:
                timer.cancel()
            return result

        return inner

    return outer