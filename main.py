from curves import Oval
import numpy as np
import matplotlib.pyplot as plt

class ConjugationPoints:

    def __init__(self):
        super(ConjugationPoints, self).__init__()

    @staticmethod
    def solve_system(points):
        count = 5
        if len(points) != count:
            raise Exception("len(points) must be {}".format(count))
        A = np.asarray([
            np.asarray([p.x**2, p.y**2, p.x*p.y, p.x, p.y])
            for p in points
        ])
        B = np.asarray([1 for i in range(count)])
        return np.linalg.solve(A, B)

    @staticmethod
    def calculate(points):
        result = []
        for i in range(2, len(points) - 2, 1):
            solution = ConjugationPoints.solve_system(points[i - 2:i + 3])
            result.append((points[i], solution[0]))
        return result

def draw(result):
    x = np.asarray([item[0].x for item in result])
    y = np.asarray([item[1] for item in result])
    plt.scatter(x, y, [2 for i in x])
    plt.show()

if __name__ == '__main__':
    oval = Oval(0.001)
    oval.draw()
    up_points = oval.points[0:int(len(oval.points)/2)]
    down_points = oval.points[int(len(oval.points) / 2) + 1:int(len(oval.points))]
    up = ConjugationPoints.calculate(up_points)
    down = ConjugationPoints.calculate(down_points)

    draw(up)
    draw(down)
