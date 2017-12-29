import math as math
import numpy as np

class circleMatrix:

    def __init__(self, diameter, depth):

        self.diameter = diameter
        self.radius = diameter/2
        self.doseData150KVP = {00.0: 1,
                               4.10: .996,
                               7.00: .992,
                               9.90: .984,
                               13.10: .971,
                               16.10: .956,
                               19.40: .938,
                               22.70: .917,
                               25.40: .894,
                               29.10: .866,
                               32.70: .837,
                               37.00: .805,
                               41.60: .768,
                               47.40: .723,
                               53.10: .676,
                               60.40: .619,
                               67.80: .566,
                               74.10: .525,
                               79.70: .489,
                               85.60: .452,
                               91.90: .414,
                               97.40: .384,
                               106.60: .341,
                               113.80: .312,
                               119.50: .291,
                               127.80: .264,
                               134.60: .240,
                               143.20: .211,
                               152.40: .184,
                               159.30: .168,
                               166.70: .152,
                               182.40: .129,
                               195.90: .111,
                               199.90: .107}
                            # x(mm): Isub0(%/100)

        self.mu = {'water': .1505, 'muscle': .1492, 'fat': .1500, 'bone': .1480}
        self.rho = {'water': 1.0, 'muscle': 1.1, 'fat': 0.9, 'bone': 2.3}
        self.gamma = self.calcGamma()

        self.intensityMatrix = np.zeros((diameter, diameter))
        self.populateMatrix()

    def calcGamma(self):
        gamma = {}
        for tissue in self.mu.keys():
            tissueGamma = self.mu[tissue] * self.rho[tissue]
            gamma.update({tissue: tissueGamma})
        return gamma

    def polarToX(self, theta):
        return self.radius * math.cos(theta)

    def polarToY(self, theta):
        return self.radius * math.sin(theta)

    def yToTheta(self, y):
        return math.asin((y / self.radius))

    def calcIntensity(self, depth, initI, tissue):
        return initI * math.exp((-1 * self.gamma[tissue]) * depth)

    def calcIntermediateSlope(self, initX, initY, finX, finY):
        deltaX = finX - initX
        deltaY = finY - initY
        slope = deltaY / deltaX
        return slope

    def calcNextSegmentLength(self, y):
        theta = self.yToTheta(y)
        halfLength = self.polarToX(theta)
        return halfLength * 2

    def populateMatrixBottomHemisphere(self):
        temp = np.zeros((self.diameter, self.diameter))
        topHemishpere = np.zeros((self.diameter, self.diameter))
        keys = []
        for key in self.doseData150KVP.keys():
            keys.append(key)

        for height in range(int(self.radius)):
            currentMinXIndex = 0
            currentMaxXIndex = 1
            currentSlope = self.calcIntermediateSlope(keys[currentMinXIndex],
                                                      self.doseData150KVP[keys[currentMinXIndex]],
                                                      keys[currentMaxXIndex],
                                                      self.doseData150KVP[keys[currentMaxXIndex]])
            maxDepth = self.calcNextSegmentLength(height)
            for depth in range(int(maxDepth)):
                while depth > keys[currentMaxXIndex] and currentMaxXIndex < len(keys):
                    currentMinXIndex = currentMaxXIndex
                    currentMaxXIndex += 1
                    currentSlope = self.calcIntermediateSlope(keys[currentMinXIndex],self.doseData150KVP[keys[currentMinXIndex]],
                                                              keys[currentMaxXIndex], self.doseData150KVP[keys[currentMaxXIndex]])
                initI = self.doseData150KVP[keys[currentMinXIndex]] + (currentSlope * (depth - keys[currentMinXIndex]))
                temp[height][depth] = initI
                topHemishpere[height][depth] = self.calcIntensity(depth, initI, 'water')
        return topHemishpere

    def populateMatrix(self):
        bottomHemisphere = self.populateMatrixBottomHemisphere()
        for height in range(int(self.radius), int(self.diameter)):
            self.intensityMatrix[height][:] = bottomHemisphere[int(height-self.radius)][:]

        for height in range(int(self.radius)):
            self.intensityMatrix[height][:] = bottomHemisphere[int(self.radius-height)][:]

        print('yay')

if __name__ == '__main__':
    cm = circleMatrix(100, 40)
    # file = open('output.txt', 'w')
    # for row in range(int(100)):
    #     file.write(str(cm.intensityMatrix[row][:]))
    # file.close()

