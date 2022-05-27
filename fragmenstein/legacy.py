"""
This file simply deals with some things that are not available in 3.7.
"""

import sys
if sys.version_info.major != 3 or sys.version_info.minor < 8:
    # singledispatchmethod allows methods to be overloaded in their first non term parameter
    # hack to enable the backport. !pip install singledispatchmethod typing_extensions

    import functools
    from singledispatchmethod import singledispatchmethod
    functools.singledispatchmethod = singledispatchmethod  # noqa this is for the older versions

    import typing
    # TypedDict is totally ace
    from typing_extensions import TypedDict, Literal, overload

    typing.TypedDict = TypedDict  # noqa
    typing.Literal = Literal  # noqa
    typing.overload = overload  # noqa
if sys.version_info.major != 3 or sys.version_info.minor < 10:
    # Unpack allows clean annotation of **kwargs!
    from typing_extensions import Unpack  # this was only "recently" added... if you get an error: upgrade package
    import typing
    typing.Unpack = Unpack
