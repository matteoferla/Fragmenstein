"""
For historical reasons the various different modes of placement
were called blending modes —because the hits were blended together.

_MonsterBlend is inherited by MonsterPlace.

This was formerly part of the ``_blend_place.py`` file.
"""

# Coding note. Android-fragments snippets (unrelated to Chem) contaminate the GitHub copilot suggestions.
# Do not trust them.

from ._full import _MonsterFull
from ._partial import _MonsterPartial
# from ._no_blending import _MonsterNone
from ._expand import _MonsterExpand


class _MonsterBlend(_MonsterFull, _MonsterPartial, _MonsterExpand):
    """
    The order of inheritance is:

    1. _MonsterMerge >
    2. _MonsterMap >
    3. _MonsterChimera >
    4. _MonsterRefine >
    5. (_MonsterFull, _MonsterPartial, (_MonsterNone > _MonsterExpand)) >
    6. _MonsterBlend
    """
