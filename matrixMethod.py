import math as math
import numpy as np

class circleMatrix:

    def __init__(self, diameter, depth):

        self.diameter = diameter
        self.radius = diameter/2
        self.innerRadius = math.fabs(self.radius - depth)
        self.zeroDepth = depth                # minimum X
        self.inverseDepth = self.radius + self.innerRadius # maximum X
        self.doseData150KVP = {00.0: 1,
                               41.0: .996,
                               70.0: .992,
                               99.0: .984,
                               131.0: .971,
                               161.0: .956,
                               194.0: .938,
                               227.0: .917,
                               254.0: .894,
                               291.0: .866,
                               327.0: .837,
                               370.0: .805,
                               416.0: .768,
                               474.0: .723,
                               531.0: .676,
                               604.0: .619,
                               678.0: .566,
                               741.0: .525,
                               797.0: .489,
                               856.0: .452,
                               919.0: .414,
                               974.0: .384, }
                            # x(mm): Isub0(%/100)

        self.mu = {'water': .1505, 'muscle': .1492, 'fat': .1500, 'bone': .1480}
        self.rho = {'water': 1.0, 'muscle': 1.1, 'fat': 0.9, 'bone': 2.3}
        self.gamma = self.calcGamma()

        self.intensityMatrix = np.zeros((diameter, diameter))

    def calcGamma(self):
        gammaDic = {}
        for tissue in self.mu.keys():
            gamma = self.mu[tissue] * self.rho[tissue]
            gammaDic.update({tissue: gamma})
        return gammaDic

    def calcIntermediateSlope(self, initX, initY, finX, finY):
        delaX = finX - initX
        deltaY = finY - initY
        slope = deltaY / delaX
        return slope

    def populateMatrix(self):
        return None

if __name__ == '__main__':
    cm = circleMatrix(100, 40)
