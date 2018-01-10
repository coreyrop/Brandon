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

    def pointDist(self, xa, ya, xb, yb):
        dist = math.sqrt((xb - xa)**2 + (yb - ya)**2)
        return dist

    def surfaceTheta(self, y):
        return self.thetaFromYR(y, self.radius)

    def calcInstantIntensity(self, depth, theta, initI, tissue):
        posCos = math.fabs(math.cos(math.radians(theta)))
        gammaFactor = math.exp((-1 * self.gamma[tissue] * depth))
        intensity = initI * posCos * gammaFactor
        return intensity

    def determineTissue(self, x, y):
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

    def determinePath(self, x, y, stepSize = 1.0, thetaOffSet = 0.0):
        path = []
        deltaX = self.polarToX(stepSize, thetaOffSet)
        deltaY = self.polarToY(stepSize, thetaOffSet)
        oldX = x
        oldY = y
        newX = x
        newY = y
        oldTissue = self.determineTissue(x, y)
        r = self.radiusFromXY(x, y)
        sectionDepth = 0
        while (r <= self.radius):
            sectionDepth += self.pointDist(oldX, oldY, newX, newY)
            newTissue = self.determineTissue(newX, newY)
            if newTissue is not oldTissue:
                path.append((oldTissue, sectionDepth))
                sectionDepth = 0
                oldTissue = newTissue
            oldX = newX
            oldY = newY
            newX += deltaX
            newY += deltaY
            r = self.radiusFromXY(newX, newY)
        path.append((oldTissue, sectionDepth))
        return path


    def calcIntensityFromPath(self, path, surfaceTheta):
        initI = 1
        intensity = None
        while(path):
            tissue, depth = path.pop()
            intensity = self.calcInstantIntensity(depth, surfaceTheta, initI, tissue)
            initI = intensity
        return Decimal(intensity) / Decimal(1)

    def populateStill(self):
        for row in range(int(self.diameter)):
            y = self.radius - row
            for col in range(int(self.diameter)):
                x = self.radius - col
                r = self.radiusFromXY(x, y)
                if r > self.radius:
                    continue
                surfaceTheta = self.surfaceTheta(y)
                path = self.determinePath(x, y)
                intensity = self.calcIntensityFromPath(path, surfaceTheta)
                self.intensityMatrixStill[row][col] = intensity

        pass
    
    def populateRotation(self):
        for row in range(int(self.diameter)):
            y = self.radius - row
            for col in range(int(self.diameter)):
                x = self.radius - col
                r = self.radiusFromXY(x, y)
                if r > self.radius:
                    continue
                for deltaTheta in range(360):
                    surfaceTheta = self.surfaceTheta(y) + deltaTheta
                    path = self.determinePath(x, y, 1, deltaTheta)
                    intensity = self.calcIntensityFromPath(path, surfaceTheta)
                    self.intensityMatrixRotated[row][col] = Decimal(self.intensityMatrixRotated[row][col]) + Decimal(intensity)
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
