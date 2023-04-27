from .base import FragmensteinParserBase
from .monster import FragmensteinParserMonster
from .victor import FragmensteinParserVictor
from .laboratory import FragmensteinParserLaboratory
from .utils import FragmensteinParserUtils


class FragmensteinParser(FragmensteinParserMonster,
                         FragmensteinParserVictor,
                         FragmensteinParserLaboratory,
                         FragmensteinParserUtils,
                         FragmensteinParserBase):
    pass
