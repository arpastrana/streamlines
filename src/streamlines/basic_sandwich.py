__name__ = "Basic Sandwich"
__author__ = "Rafael Pastrana"
__version__ = "0.0.2"
__date__ = "2018.08.24"


"""
Reinforcement design of concrete shells acc. to the FIB Bulletin 45 - 2008.

It is based on the basic 3-layer sandwich model developed by Lourenco and
it accounts for both membrane and flexural forces.

Given a coordinate system XY, this module proposes reinforcement in 4 layers:
    - 2 layers are hosted at the Top Layer of the Shell
    - 2 layers are hosted at the Bottom Layer of the Shell

The reinforcement is assumed to be parallel to the XY directions of the
coordinate system. In other words, they comprise an orthogonal arrangement.

Known limitations of the basic 3-layer sandwich model are:
    - The shell core doesn't contribute to transfer membrane forces
    - The same lever arm is assumed for both the Top and Bottom layers of rebar
    - The angle of concrete crack at the core is taken as 45 degrees by default
    - Biaxial compression is disregarded.

Such limitations are addresed by the advanced 3-layer sandwich model developed
by Lourenco himself. The advanced model is explained on J. Blauauwendraad's
book "Plates and FEM: Surprises and Pitfalls" from 2010.
"""


import math


class ReinforcedShell(object):

    def __init__(self, id, matConcrete, matRebar,
                 height=0.1, cover=0.02, theta=math.radians(45),
                 ):

        # 0. Meta Attributes
        self.id = id

        # 1. Geometric Properties
        self.height = height
        self.cover = cover
        self.tbot = self.cover
        self.tsup = self.cover
        self.dv = self.height - self.cover
        self.theta = theta
        self.thetaDeg = math.degrees(self.theta)
        self.ktop = 1
        self.kbot = 1
        self.kcen = 1

        # 2. Material Properties
        self.concrete = matConcrete
        self.rebar = matRebar
        self.fcd = self.concrete.fcd
        self.fctd = self.concrete.fctd
        self.fyd = self.rebar.fyd

        # 3. Acting Forces
        self.nx = 0.0     # membrane force on the X-Axis in kN/m
        self.ny = 0.0     # membrane force on the Y-axis in kN/m
        self.nxy = 0.0    # membrane shear force in kN/m

        self.mx = 0.0    # flexural bending moment around the X-axis
        self.my = 0.0    # flexural bending moment around the Y-axis
        self.mxy = 0.0   # torsional bending moment

        self.vx = 0.0    # transversal shear force along the X-axis
        self.vy = 0.0    # transversal shear force along  the Y-axis

        # 4. Crush Flags
        self.coreCrush = False
        self.tsupCrush = False
        self.tinfCrush = False

        # 5. Resistance Forces
        self.nxt = 0.0      # force at the top along the X-Axis in kN/m
        self.nxb = 0.0      # force at the bottom along the X-Axis in kN/m
        self.nxc = 0.0      # force at the center along the X-Axis in kN/m

        self.nyt = 0.0      # force at the top along the Y-Axis in kN/m
        self.nyb = 0.0      # force at the bottom along the Y-Axis in kN/m
        self.nyc = 0.0      # force at the center along the Y-Axis in kN/m

        self.nxyt = 0.0     # transverse shear force at the top in kN/m
        self.nxyb = 0.0     # transverse shear force at the bottom in kN/m
        self.nxyc = 0.0     # transverse shear force at the center in kN/m

        self.tcred = 0.30 * self.concrete.fctd * (1.6 - self.dv)    # MPa
        # ENV 1992-1-1
        # assuming not high tensile forces
        # how to know how high is a high tensile force?
        # it disregards positive effect of long. reinforcement
        # it disregards in-plane compressive stress
        # self.tc = 0.25 * self.rebar.fctd * (1.6 - self.dv) *
        # (1.2 + 40 * self.rhol) + 0.15 * self.sigmacp

        # 6. Reinforcement design
        self.asxt = -1  # area of rebar at the top layer along the X-axis
        self.asxb = -1  # area of rebar at the bottom layer along the X-axis
        self.asxc = -1
        self.asyt = -1  # area of rebar at the top layer along the Y-axis
        self.asyb = -1  # area of rebar at the bottom layer along the Y-axis
        self.asyc = -1
        self.ast = -1   # area of transverse reinforcement

    def checkCoreCrush(self):
        if self.vo * 0.001 / self.dv > self.tcred:
            self.coreCrush = True   # core is craked

    def checkCoverCrush(self):
        a = 0.001
        # check cover crush at top layer
        fct = ((self.nsxt + self.nsyt) - (self.nxt + self.nyt)) * a / self.tsup
        if fct > self.concrete.fcd:
            self.tsupCrush = True
        # check cover crush at bottom layer
        fcb = ((self.nsxb + self.nsyb) - (self.nxb + self.nyb)) * a / self.tbot
        if fcb > self.concrete.fcd:
            self.tinfCrush = True

    def applyForces(self, nx, ny, nxy, mx, my, mxy, vx, vy):
        self.nx = nx     # membrane force on the X-Axis in kN/m
        self.ny = ny     # membrane force on the Y-axis in kN/m
        self.nxy = nxy    # membrane shear force in kN/m

        self.mx = mx   # flexural bending moment around the X-axis
        self.my = my    # flexural bending moment around the Y-axis
        self.mxy = mxy   # torsional bending moment

        self.vx = vx   # transversal shear force along the X-axis
        self.vy = vy    # transversal shear force along  the Y-axis

        # transversal principal shear force
        self.vo = math.sqrt(self.vx * self.vx + self.vy * self.vy)
        self.phi = math.atan(self.vx / self.vy)
        self.phiDeg = math.degrees(self.phi)

    def getMembraneForces(self):

        self.checkCoreCrush()

        # membrane forces - top and bottom - kN/m
        nxt = (-1.0 * self.mx / self.dv) + self.nx / 2
        nxb = (1.0 * self.mx / self.dv) + self.nx / 2
        nyt = (-1.0 * self.my / self.dv) + self.ny / 2
        nyb = (1.0 * self.my / self.dv) + self.ny / 2

        # membrane shear forces - top and bottom - kN/m
        nxyt = (-1.0 * self.mxy / self.dv) + self.nxy / 2
        nxyb = (1.0 * self.mxy / self.dv) + self.nxy / 2

        # transverse shear forces
        nvx = 0.0
        nvy = 0.0
        nvxy = 0.0

        if self.coreCrush is True:
            # transverse shear forces - if core is cracked, replace - kN/m
            nvx = self.vx * self.vx / (2 * self.vo * math.tan(self.theta))
            nvy = self.vy * self.vy / (2 * self.vo * math.tan(self.theta))
            nvxy = self.vx * self.vy / (2 * self.vo * math.tan(self.theta))

        self.nxt = nxt + nvx
        self.nxb = nxb + nvx

        self.nyt = nyt + nvy
        self.nyb = nyb + nvy

        self.nxyt = nxyt + nvxy
        self.nxyb = nxyb + nvxy

    def get_membrane_forces(self):

        self.checkCoreCrush()

        # membrane forces - top and bottom - kN/m
        # nxt = (-1.0 * self.mx / self.dv) + self.nx / 2
        # nxb = (1.0 * self.mx / self.dv) + self.nx / 2
        # nyt = (-1.0 * self.my / self.dv) + self.ny / 2
        # nyb = (1.0 * self.my / self.dv) + self.ny / 2

        nxt = (-1.0 * self.mx / self.dv)
        nxb = (1.0 * self.mx / self.dv)
        nyt = (-1.0 * self.my / self.dv)
        nyb = (1.0 * self.my / self.dv)

        # membrane axial forces - center - kN/m
        nxc = self.nx
        nyc = self.ny

        # membrane shear forces - top and bottom - kN/m
        # nxyt = (-1.0 * self.mxy / self.dv) + self.nxy / 2
        # nxyb = (1.0 * self.mxy / self.dv) + self.nxy / 2

        nxyt = (-1.0 * self.mxy / self.dv)
        nxyb = (1.0 * self.mxy / self.dv)

        # shear forces center
        nxyc = self.nxy

        # transverse shear forces
        # nvx = 0.0
        # nvy = 0.0
        # nvxy = 0.0

        # if self.coreCrush is True:
        #     # transverse shear forces - if core is cracked, replace - kN/m
        #     nvx = self.vx * self.vx / (2 * self.vo * math.tan(self.theta))
        #     nvy = self.vy * self.vy / (2 * self.vo * math.tan(self.theta))
        #     nvxy = self.vx * self.vy / (2 * self.vo * math.tan(self.theta))

        self.nxt = nxt
        self.nxb = nxb
        self.nxc = nxc

        self.nyt = nyt
        self.nyb = nyb
        self.nyc = nyc

        self.nxyt = nxyt
        self.nxyb = nxyb
        self.nxyc = nxyc

    def getReinforcementForces(self):
        # forces at the bottom
        self.nsxb = (self.nxb + self.kbot * math.fabs(self.nxyb))
        self.nsyb = (self.nyb + math.pow(self.kbot, -1) * math.fabs(self.nxyb))

        # forces at the top
        self.nsxt = self.nxt + self.ktop * math.fabs(self.nxyt)
        self.nsyt = self.nyt + math.pow(self.ktop, -1) * math.fabs(self.nxyt)

        # transverse shear force
        self.nst = self.vo * math.tan(self.theta)

    def get_reinforcement_forces(self):
        # forces at the bottom
        self.nsxb = self.nxb + self.kbot * math.fabs(self.nxyb)
        self.nsyb = self.nyb + math.pow(self.kbot, -1) * math.fabs(self.nxyb)

        # forces at the top
        self.nsxt = self.nxt + self.ktop * math.fabs(self.nxyt)
        self.nsyt = self.nyt + math.pow(self.ktop, -1) * math.fabs(self.nxyt)

        # force at the middle
        self.nsxc = self.nxc + self.kcen * math.fabs(self.nxyc)
        self.nsyc = self.nyc + math.pow(self.kcen, -1) * math.fabs(self.nxyc)

        # transverse shear force
        self.nst = self.vo * math.tan(self.theta)

    def getReinforcementBottom(self):
        a = 10.0
        self.asxb = self.nsxb * a / self.rebar.fyd
        self.asyb = self.nsyb * a / self.rebar.fyd

        if self.asxb < 0.0:
            self.asxb = 0.0
        if self.asyb < 0.0:
            self.asyb = 0.0

    def getReinforcementTop(self):
        a = 10.0
        self.asxt = self.nsxt * a / self.rebar.fyd
        self.asyt = self.nsyt * a / self.rebar.fyd

        if self.asxt < 0.0:
            self.asxt = 0.0
        if self.asyt < 0.0:
            self.asyt = 0.0

    def getReinforcementCenter(self):
        a = 10.0
        self.asxc = self.nsxc * a / self.rebar.fyd
        self.asyc = self.nsyc * a / self.rebar.fyd

        if self.asxc < 0.0:
            self.asxc = 0.0
        if self.asyc < 0.0:
            self.asyc = 0.0

    def getReinforcementTransverse(self):
        a = 10.0
        if self.coreCrush is True:
            self.ast = self.nst * a / self.rebar.fyd
        else:
            self.ast = 0.0

    def designReinforcement(self):
        # 1. Check if Covers are Crushed. If they are, don't perform design.
        # Steel Areas asxt, asyt, asxb, and asyb would be kept as -1 cm2/m
        self.checkCoverCrush()

        # 2. Perform Reinforcement Design
        if self.tsupCrush is False and self.tinfCrush is False:
            # 3. Get Reinforcement at the Top Cover
            self.getReinforcementTop()
            # 4. Get Reinforcement at the Bottom Cover
            self.getReinforcementBottom()
            # 5. Compute Shear Reinforcement
            self.getReinforcementTransverse()
        # else:
        #    print("Covers have crushed. Please increase shell height.")

    def design_reinforcement(self):
        # 1. Check if Covers are Crushed. If they are, don't perform design.
        # Steel Areas asxt, asyt, asxb, and asyb would be kept as -1 cm2/m
        self.checkCoverCrush()

        # 2. Perform Reinforcement Design
        if self.tsupCrush is False and self.tinfCrush is False:
            # 3. Get Reinforcement at the Top Cover
            self.getReinforcementTop()
            # 4. Get Reinforcement at the Bottom Cover
            self.getReinforcementBottom()
            # 5. Get Reinforcement at the Center
            self.getReinforcementCenter()
            # 5. Compute Shear Reinforcement
            self.getReinforcementTransverse()
        # else:
        #    print("Covers have crushed. Please increase shell height.")

    def printSummary(self):
        print(
            """
            ###################################
            Shell ID: {id}
            -----------------------------------
            fcd = {fcd} MPa
            fctd = {fctd} MPa
            fyd = {fyd} MPa
            tcred = {tcred} MPa
            -----------------------------------
            h = {height} m
            tinf = {tbot} m
            tsup = {tsup} m
            dv = {dv} m
            -----------------------------------
            nx = {nx} kN/m
            ny = {ny} kN/m
            nxy = {nxy} kN/m
            mx = {mx} kNm/m
            my = {my} kNm/m
            mxy = {mxy} kNm/m
            vx = {vx} kN/m
            vy = {vy} kN/m
            vo = {vo} kN/m
            theta = {thetaDeg} deg
            phi = {phiDeg} deg
            -----------------------------------
            coreCrush = {coreCrush}
            tsupCrush = {tsupCrush}
            tinfCrush = {tinfCrush}
            -----------------------------------
            nxt = {nxt} kN/m
            nyt = {nyt} kN/m
            nxyt = {nxyt} kN/m
            nxb = {nxb} kN/m
            nyb = {nyb} kN/m
            nxyb = {nxyb} kN/m
            -----------------------------------
            nsxt = {nsxt} kN/m
            nsyt = {nsyt} kN/m
            nsxb = {nsxb} kN/m
            nsyb = {nsyb} kN/m
            nst = {nst} kN/m
            -----------------------------------
            asxt = {asxt} cm2/m
            asyt = {asyt} cm2/m
            asxb = {asxb} cm2/m
            asyb = {asyb} cm2/m
            ast = {ast} cm2/m
            ###################################
            """.format(**vars(self))
        )


if __name__ == "__main__":
    print('hello')
