
import math as math

class new:

    def __init__(self, diameter, depth):

        self.diameter = diameter
        self.radius = diameter/2
        self.innerRadius = math.fabs(self.radius - depth)
        self.zeroDepth = depth                # minimum X
        self.inverseDepth = self.radius + self.innerRadius # maximum X
        self.doseData150KVP = {0.00: 1,
                               .410: .996,
                               .700: .992,
                               .990: .984,
                               1.31: .971,
                               1.61: .956,
                               1.94: .938,
                               2.27: .917,
                               2.54: .894,
                               2.91: .866,
                               3.27: .837,
                               3.70: .805,
                               4.16: .768,
                               4.74: .723,
                               5.31: .676,
                               6.04: .619,
                               6.78: .566,
                               7.41: .525,
                               7.97: .489,
                               8.56: .452,
                               9.19: .414,
                               9.74: .384, }
                            # x(cm): Isub0(%/100)

        self.mu = {'water': .1505, 'muscle': .1492, 'fat': .1500, 'bone': .1480}
        self.rho = {'water': 1.0, 'muscle': 1.1, 'fat': 0.9, 'bone': 2.3}
        self.gamma = self.calcGamma()

    def calcGamma(self):
        gammaDic = {}
        for tissue in self.mu.keys():
            gamma = self.mu[tissue] * self.rho[tissue]
            gammaDic.update({tissue: gamma})
        return gammaDic

    def summationA(self):
        sumTotal = 0
        startXIntensity = None
        endXIntensity = None
        for x,I in self.doseData150KVP.items():
            if (x >= self.zeroDepth and x <= self.inverseDepth):
                currentI = I * math.exp((-1 * self.gamma['water']) * x)
                if startXIntensity is None:
                    startXIntensity = currentI
                endXIntensity = currentI
                sumTotal += currentI
        sumTotal = ((sumTotal * 2) - startXIntensity) - endXIntensity
        return sumTotal





class oneFiftyKVP:

    def __init__(self, x, y):

        self.intensity = None
        self.totalIntensity = None
        self.diameter = 10  #(cm)
        self.x = x
        self.y = y
        self.r = self.cartisianToPolar()
        self.doseData150KVP = {0.00 : 1,
                               .410 : .996,
                               .700 : .992,
                               .990 : .984,
                               1.31 : .971,
                               1.61 : .956,
                               1.94 : .938,
                               2.27 : .917,
                               2.54 : .894,
                               2.91 : .866,
                               3.27 : .837,
                               3.70 : .805,
                               4.16 : .768,
                               4.74 : .723,
                               5.31 : .676,
                               6.04 : .619,
                               6.78 : .566,
                               7.41 : .525,
                               7.97 : .489,
                               8.56 : .452,
                               9.19 : .414,
                               9.74 : .384,}
                            # x(cm) : Isub0(%/100)

        self.mu = {'water' : .1505, 'muscle' : .1492, 'fat' : .1500, 'bone' : .1480}
        self.rho = {'water' : 1.0, 'muscle' : 1.1, 'fat' : 0.9, 'bone' : 2.3}
        self.gamma = self.calcGamma()

    def polarToX(self, r, theta):
        x = r * math.cos(theta)
        return x

    def polarToY(self, r, theta):
        y = r * math.sin(theta)
        return y

    def cartisianToPolar(self):
        r = math.sqrt((self.x * self.x) + (self.y * self.y))
        return r

    def calcGamma(self):
        gammaDic = {}
        for tissue in self.mu.keys():
            gamma = self.mu[tissue] * self.rho[tissue]
            gammaDic.update({tissue : gamma})
        return gammaDic

    def calcEndX(self):
        gaps = self.diameter - (2 * self.r)
        singleGap = gaps / 2
        endX = 2 * self.r + singleGap
        return endX

    def linearCalc(self):
        endX = self.calcEndX()
        return self.sumAllIForX(self.x, endX)


    def sumAllIForX(self, startX, endX):
        totalIntensity = 0
        firstIntensity = None
        pairs = self.doseData150KVP.items()
        for x,Io in pairs:
            if (x >= startX and x <= endX) or (x <= startX and x >= endX):
                if firstIntensity is None:
                    firstIntensity = Io
                totalIntensity += Io * math.exp((-1 * self.gamma['water']) * x)
        return (totalIntensity * 2) - firstIntensity

if __name__ == '__main__':
    c = oneFiftyKVP(6,2)
    print( 'x: ' + str(c.x) + ' total Intensity: ' + str(c.linearCalc()))
    k = new(10, 4)
    print('zeroDepth: ' + str(k.zeroDepth) + ' total Intensity: ' + str(k.summationA()))