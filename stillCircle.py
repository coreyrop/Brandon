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
        self.defaultTissue = 'water'
        self.boneRegion = { 'tissue': 'bone', 'xCenter': -30, 'yCenter': 20, 'radius': 15}

        self.intensityMatrixStill = np.zeros((self.diameter, self.diameter))
        self.intensityMatrixRotated = np.zeros((self.diameter, self.diameter))
        self.populateStill()
        self.populateRotation()

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

    def calcIntermediateSlope(self, initX, initY, finX, finY):
        deltaX = finX - initX
        deltaY = finY - initY
        slope = deltaY / deltaX
        return slope

    def calcIntensity(self, depth, theta, initI, tissue):
        posCos = math.fabs(math.cos(math.radians(theta)))
        gammaFactor = math.exp((-1 * self.gamma[tissue]))
        intensity = initI * posCos * gammaFactor * depth
        return intensity

    def defineTissue(self, x, y):
        xOffOrigin = self.boneRegion['xCenter']
        yOffOrigin = self.boneRegion['yCenter']
        regionRadius = self.boneRegion['radius']
        if y > yOffOrigin + regionRadius or y < yOffOrigin - regionRadius:
            return self.defaultTissue
        relativeY = y - yOffOrigin
        thetaInRegion = self.thetaFromYR(relativeY, regionRadius)
        deltaXInRegion = self.polarToX(regionRadius, thetaInRegion)
        minXInRegion = xOffOrigin - deltaXInRegion
        maxXInRegion = xOffOrigin + deltaXInRegion
        if x <= maxXInRegion and x >= minXInRegion:
            return self.boneRegion['tissue']
        else:
            return self.defaultTissue

    def populateStill(self):
        depths = []
        for key in self.doseData150KVP.keys():
            depths.append(key)

        for row in range(int(self.diameter)):
            currentMinIndex = 0
            currentMaxIndex = 1
            y = self.radius - row
            for col in range(int(self.diameter)):
                x = self.radius - col
                resultantRadius = self.radiusFromXY(x,y)
                if resultantRadius > self.radius:
                    continue
                theta = self.thetaFromXY(x, y)
                surfaceTheta = self.thetaFromYR(y, self.radius)
                surfaceX = self.polarToX(self.radius, surfaceTheta)
                depth = surfaceX - x
                while (depth > depths[currentMaxIndex] and currentMaxIndex < len(depths)):
                    currentMinIndex = currentMaxIndex
                    currentMaxIndex += 1

                slope = self.calcIntermediateSlope(depths[currentMinIndex], self.doseData150KVP[depths[currentMinIndex]],
                                                   depths[currentMaxIndex], self.doseData150KVP[depths[currentMaxIndex]])

                initI = self.doseData150KVP[depths[currentMinIndex]] + (slope * (depth - depths[currentMinIndex]))
                tissue = self.defineTissue(x, y)
                if tissue is 'bone':
                    print('x: ' + str(x) + ' y: ' + str(y))
                i = Decimal(self.calcIntensity(depth, surfaceTheta, initI, tissue)) / Decimal(1)
                self.intensityMatrixStill[row][col] = i

        pass

    def populateRotation(self):
        depths = []
        for key in self.doseData150KVP.keys():
            depths.append(key)

        for deltaTheta in range(360):
            for row in range(int(self.diameter)):
                currentMinIndex = 0
                currentMaxIndex = 1
                y = self.radius - row
                for col in range(int(self.diameter)):
                    x = self.radius - col
                    resultantRadius = self.radiusFromXY(x,y)
                    if resultantRadius > self.radius:
                        continue
                    theta = self.thetaFromXY(x, y)
                    newTheta = theta + deltaTheta
                    newY = self.polarToY(resultantRadius, newTheta)
                    newX = self.polarToX(resultantRadius, newTheta)
                    surfaceTheta = self.thetaFromYR(newY, self.radius)
                    surfaceX = self.polarToX(self.radius, surfaceTheta)
                    depth = surfaceX - newX
                    while (depth > depths[currentMaxIndex] and currentMaxIndex < len(depths)):
                        currentMinIndex = currentMaxIndex
                        currentMaxIndex += 1

                    slope = self.calcIntermediateSlope(depths[currentMinIndex], self.doseData150KVP[depths[currentMinIndex]],
                                                       depths[currentMaxIndex], self.doseData150KVP[depths[currentMaxIndex]])

                    initI = self.doseData150KVP[depths[currentMinIndex]] + (slope * (depth - depths[currentMinIndex]))
                    tissue = self.defineTissue(x, y)

                    i = Decimal(self.calcIntensity(depth, surfaceTheta, initI, tissue)) / Decimal(1)
                    self.intensityMatrixRotated[row][col] = Decimal(self.intensityMatrixRotated[row][col]) + Decimal(i)

        pass

if __name__ == '__main__':
    circMatrx = circleMatrix(10, 1, 4)
    print('yay')

    # 2D plot of intensity vs. depth across diameter
    ####################################################################################################################
    # fig, ax = plt.subplots()
    #
    # depth = []
    # intensity = circMatrx.intensityMatrixStill[int(circMatrx.radius)][:]
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

    # 3D plot of intensity vs. depth [still model]
    ####################################################################################################################
    depth = [x for x in range(circMatrx.diameter)]
    height = [y for y in range(circMatrx.diameter)]
    intensity = circMatrx.intensityMatrixStill

    depth, height = np.meshgrid(depth, height)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_title('Still model intensity plot')
    ax.set_xlabel('depth (mm)')
    ax.set_ylabel('height (mm)')
    ax.set_zlabel('intensity (% / 100)')
    ax.plot_surface(depth, height, intensity, rstride=1, cstride=1, cmap=cm.viridis)


    plt.show()
    ####################################################################################################################

    # 3D plot of intensity vs. depth [rotated model]
    ####################################################################################################################
    depth = [x for x in range(circMatrx.diameter)]
    height = [y for y in range(circMatrx.diameter)]
    intensity = circMatrx.intensityMatrixRotated

    depth, height = np.meshgrid(depth, height)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_title('Rotation model intensity plot')
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
