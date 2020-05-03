from enum import Enum


class EtatVol(Enum):
    SOL = 0  # FIXME : SOL != RAMPE ? (Pour l'instant il sert a rien cet état, c'était dans l'optique où on fait une simu embarquée en temps réel)
    RAMPE = 1
    VOL = 2
    SAUVEGARDE = 3
    ATTERI = 4
