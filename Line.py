import math

class Line():
    """Represents a 2D line that one end of a Beam will slide along
    slope-intercept form: ax+by+c=0
    """
    zeroThreshold = .00001
    
    
    def __init__(self, a, b, c):
        self.a = a+0.0
        self.b = b+0.0
        self.c = c+0.0

    def distanceToPoint(self, x, y):
        return self.a*x+self.b*y+self.c/math.sqrt(self.a**2+self.b**2)

    def closestPoint(self, x, y):
        xLine = (self.b(self.b*x-self.a*y)-self.a*self.c)/(self.a**2+self.b**2)
        yLine = (self.a(-self.b*x+self.a*y)-self.b*self.c)/(self.a**2+self.b**2)

    def slope(self):
        return -self.a/self.b

    def solveForY(self, x):
        return (-self.c-self.a*x)/self.b

    def solveForX(self, y):
        return (-self.c-self.b*y)/self.a