from .base import FragmensteinParserBase
from .monster import FragmensteinParserMonster
from .victor import FragmensteinParserVictor
from .laboratory import FragmensteinParserLaboratory
from .pipeline import FragmensteinParserPipeline
from .utils import FragmensteinParserUtils


class FragmensteinParser(FragmensteinParserMonster,
                         FragmensteinParserVictor,
                         FragmensteinParserLaboratory,
                         FragmensteinParserPipeline,
                         FragmensteinParserUtils,
                         FragmensteinParserBase):
    pass
