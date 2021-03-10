import functools
import threading

'''
https://stackoverflow.com/questions/308999/what-does-functools-wraps-do
'''
def timeout(timeout, raise_exc=False):
    """
    raise_exc - if exception should be raised on timeout
                or exception inside decorated func.
                Otherwise None will be returned.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            res = None
            exc = None
            def _run():
                nonlocal res
                nonlocal exc
                try:
                    res = func(*args, **kwargs)
                except Exception as e:
                    exc = e
            t = threading.Thread(target=_run)
            t.daemon = True
            t.start()
            t.join(timeout=timeout)
            if raise_exc and t.is_alive():
                raise TimeoutError()
            elif raise_exc and (exc is not None):
                raise exc
            else:
                return res
        return wrapper
    return decorator

def test():
    import time
    @timeout(2)
    def dummy(x):
        time.sleep(4)
        return 1

    dummy(3)


if __name__ == "__main__":
    test()