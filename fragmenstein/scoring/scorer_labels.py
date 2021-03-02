SCORE_NAME_TEMPLATE = "%s_score"
FRAGMENTS_ID = "fragments"
MOL_NAME_ID = "name"


def checkIfNameIsScore(name):
    return name.endswith(SCORE_NAME_TEMPLATE % "")

def removeScoreTag(name):
    return name.replace(SCORE_NAME_TEMPLATE % "", "")