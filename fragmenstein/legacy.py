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
if sys.version_info.major != 3 or sys.version_info.minor < 11:
    # Unpack allows clean annotation of **kwargs!
    from typing_extensions import Unpack, NotRequired
    # this was only "recently" added... if you get an error: upgrade package typing_extensions
    import typing
    typing.Unpack = Unpack
    typing.NotRequired = NotRequired


# NGLView has been broken for a while. This is a hack to ignore it.
# the issue stems from widgets 7 to 8.
try:
    import nglview
except Exception as err:
    from unittest.mock import MagicMock
    import sys

    sys.modules['nglview'] = MagicMock()
    sys.modules['nglview.component'] = MagicMock()
