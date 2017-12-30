import math as math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from decimal import getcontext, Decimal

class circleMatrix:

    def __init__(self, diameter, precisionDepth, precisionIntensity):

        getcontext().prec = precisionIntensity
        self.diameter = diameter * (10 * precisionDepth)
        self.radius = self.diameter/2
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

        self.intensityMatrix = np.zeros((self.diameter, self.diameter))
        self.populateMatrix()
        self.rotateMatrixSummation()

    def calcGamma(self):
        gamma = {}
        for tissue in self.mu.keys():
            tissueGamma = self.mu[tissue] * self.rho[tissue]
            gamma.update({tissue: tissueGamma})
        return gamma

    def polarToX(self, r, theta):
        return r * math.cos(math.radians(theta))

    def polarToY(self, r, theta):
        return r * math.sin(math.radians(theta))

    def thetaFromYR(self, y, r):
        if r == 0:
            return 0
        return math.degrees(math.asin((y / r)))

    def thetaFromXR(self, x, r):
        if r == 0:
            return 0
        return math.degrees(math.acos((x / r)))

    def thetaFromXY(self, x, y):
        if x == 0:
            return 0
        return math.degrees(math.atan2(y , x))

    def radiusFromXY(self, x, y):
        return math.sqrt((x**2) + (y**2))

    def calcIntensity(self, depth, initI, tissue):
        return  Decimal(initI * math.exp((-1 * self.gamma[tissue]) * depth)) / Decimal(1)

    def calcIntermediateSlope(self, initX, initY, finX, finY):
        deltaX = finX - initX
        deltaY = finY - initY
        slope = deltaY / deltaX
        return slope

    def calcNextSegmentLength(self, height):
        if height > self.radius:
            heightPastRadius = height - self.radius
            height = self.radius - heightPastRadius
        theta = self.thetaFromYR(height, self.radius)
        halfLength = self.polarToX(self.radius, theta)
        return halfLength * 2

    def populateMatrixBottomHemisphere(self):
        temp = np.zeros((self.diameter, self.diameter))
        bottomHemishpere = np.zeros((self.diameter, self.diameter))
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
                bottomHemishpere[height][depth] = self.calcIntensity(depth, initI, 'water')
        return bottomHemishpere

    def populateMatrix(self):
        bottomHemisphere = self.populateMatrixBottomHemisphere()
        for height in range(int(self.radius), int(self.diameter)):
            self.intensityMatrix[height][:] = bottomHemisphere[int(height-self.radius)][:]

        for height in range(int(self.radius)):
            self.intensityMatrix[height][:] = bottomHemisphere[int(self.radius-height)][:]

        print('yay')

    def rotateMatrixSummation(self):
        rotationMatrix = self.intensityMatrix.copy()
        for deltaTheta in range(1, 360):

            for height in range(int(self.diameter)):
                y = self.radius - height
                maxDepth = self.calcNextSegmentLength(math.fabs(y))
                halfMaxDepth = maxDepth / 2
                startingX = (-1 * self.radius) + (self.radius - halfMaxDepth)
                for depth in range(int(maxDepth)):
                    x = startingX + depth
                    r = self.radiusFromXY(x, y)
                    thetaxy = self.thetaFromXY(x, y)
                    calcy = self.polarToY(r, thetaxy)
                    calcx = self.polarToX(r, thetaxy)
                    newTheta = thetaxy + deltaTheta
                    newY = self.polarToY(r, newTheta)
                    newX = self.polarToX(r, newTheta)
                    newR = int(self.radiusFromXY(newX, newY))
                    newMaxDepth = self.calcNextSegmentLength(math.fabs(newY))
                    newHalfMaxDepth = newMaxDepth / 2
                    newStartingX = (-1 * self.radius) + (self.radius - newHalfMaxDepth)
                    deltaX = newX - (newStartingX + depth)
                    deltaY = newY + y
                    newDepth = int(depth + deltaX)
                    newHeight = int(self.radius - newY)
                    if newDepth >= self.diameter:
                        newDepth = self.diameter-1
                    if newHeight >= self.diameter:
                        newHeight = self.diameter-1
                    self.intensityMatrix[height][depth] = Decimal(self.intensityMatrix[height][depth]) + Decimal(rotationMatrix[newHeight][newDepth]) # with 3 sigfig precision
                    # self.intensityMatrix[height][depth] += rotationMatrix[newHeight][newDepth]                                                        # with full precision
        pass



if __name__ == '__main__':
    circMatrx = circleMatrix(10, 1, 3)
    print('yay')

    # 2D plot of intensity vs. depth across diameter
    ####################################################################################################################
    # fig, ax = plt.subplots()
    #
    # depth = []
    # intensity = circMatrx.intensityMatrix[int(circMatrx.radius)][:]
    # for depthVal in range(circMatrx.diameter):
    #     depth.append(depthVal)
    #
    # ax.plot(depth, intensity, 'k--', label='intensity as a function of depth')
    # ax.set_xlabel('depth (mm)')
    # ax.set_ylabel('intensity (% / 100)')
    # legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
    # legend.get_frame().set_facecolor('#00FFCC')
    # plt.show()
    ####################################################################################################################

    # 3D plot of intensity vs. depth
    ####################################################################################################################
    depth = [x for x in range(circMatrx.diameter)]
    height = [y for y in range(circMatrx.diameter)]
    intensity = circMatrx.intensityMatrix

    depth, height = np.meshgrid(depth, height)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('depth (mm)')
    ax.set_ylabel('height (mm)')
    ax.set_zlabel('intensity (% / 100)')
    ax.plot_surface(depth, height, intensity, rstride=1, cstride=1, cmap=cm.viridis)


    plt.show()
    ####################################################################################################################

    # 3D plot example
    ####################################################################################################################
    # X = np.arange(-5, 5, 0.25)
    # Y = np.arange(-5, 5, 0.25)
    # X, Y = np.meshgrid(X, Y)
    # R = np.sqrt(X ** 2 + Y ** 2)
    # Z = np.sin(R)
    #
    # fig = plt.figure()
    # ax = Axes3D(fig)
    # ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.viridis)
    #
    # plt.show()

