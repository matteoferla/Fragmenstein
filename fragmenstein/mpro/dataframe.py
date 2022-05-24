# The chemists at Postera were in charge of the data management
# this should probably called something different

import io

import pandas as pd
import requests


def fetch_postera() -> pd.DataFrame:
    """
    Reads the submission file off Github.
    For a local version, just ``postera = read_postera(file)``.
    :return:
    """
    url = "https://raw.githubusercontent.com/postera-ai/" + \
          "COVID_moonshot_submissions/master/covid_submissions_all_info.csv"
    s = requests.get(url).content
    postera = pd.read_csv(io.StringIO(s.decode('utf-8')))
    _add_category(postera)
    return postera


def read_postera(filename: str) -> pd.DataFrame:
    postera = pd.read_csv(filename)
    _add_category(postera)
    return postera


# not exported
def _get_category(row: pd.Series) -> str:
    """
    Postera table has categories as True/False. But it is unlikely that there are multiple.
    Turns out these categories are **not** user submitted.
    However, for consistency with other analysis by other people these are used.

    :param postera: pandas table modified in place
    :return:
    """
    for category in ('Acrylamide', 'Chloroacetamide', 'Vinylsulfonamide', 'Nitrile'):
        if row[category] in ('True', 'true', True):
            return category
    else:
        return 'non-covalent'


def _add_category(postera: pd.DataFrame) -> None:
    postera['category'] = postera.apply(_get_category, axis=1)
