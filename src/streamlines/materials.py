__name__ = "Materials"
__author__ = "Rafael Pastrana"
__version__ = "0.0.2"
__creation__ = "2018.11.12"
__date__ = "2018.11.12"


import math


class Material(object):
    def __init__(self, name="Blank Material"):
        self.name = name


class Concrete(Material):
    def __init__(self, fck=40):
        Material.__init__(self, name="Concrete")
        self.fck = fck    # fck: char. cylinder compression strength - MPa
        self.gamma = 1.5
        self.fcm = self.fck + 8 # fcm: mean compressive strength - MPa
        self.fcd = 1.0 * self.fck / self.gamma # design compression strength - MPa

        if self.fck <= 50:
            self.fctm = 0.3 * math.pow(self.fck, 0.6667)  # mean tensile strength - MPa
        elif self.fck > 50:
            temp = 1 + (self.fcm / 10)
            self.fctm = 2.12 * math.log1p(temp) # mean tensile strength - MPa
        self.fctk5 = 0.7 * self.fctm   # 5% fractile tensile strenght - MPa
        self.fctd = 1.0 * self.fctk5 / self.gamma

    def printSummary(self):
        print(
        """
        ###################################
        Material Name: {name}
        -----------------------------------
        fck = {fck} MPa
        fcm = {fcm} MPa
        fcd = {fcd} MPa
        fctm = {fctm} MPa
        fctk5 = {fctk5} MPa
        fctd = {fctd} MPa
        -----------------------------------
        """.format(**vars(self)))


class Rebar(Material):
    def __init__(self,fy=500):
        self.fy = fy    # units: MPa
        self.gamma = 1.15
        self.fyd = 1.0 * self.fy / self.gamma
