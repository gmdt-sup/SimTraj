from enum import Enum


class EtatVol(Enum):
    SOL = 0 # FIXME : SOL != RAMPE ?
    RAMPE = 1
    VOL = 2
    SAUVEGARDE = 3
    ATTERI = 4
