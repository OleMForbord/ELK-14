import sys
sys.path.append(".")
from math import *


def tij(
    i, j, teta, g, b
):  # Creates the tij needed in the derivatives. i,j=bus nr i and j.teta=the angles g=conductance matrix, b=susceptance matrix
    return g[i - 1][j - 1] * cos(teta[i - 1] - teta[j - 1]) + b[i - 1][j - 1] * sin(
        teta[i - 1] - teta[j - 1]
    )


def uij(i, j, teta, g, b):  # Creeates the uij needed in the derivatives
    return g[i - 1][j - 1] * sin(teta[i - 1] - teta[j - 1]) - b[i - 1][j - 1] * cos(
        teta[i - 1] - teta[j - 1]
    )
