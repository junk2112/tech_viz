from curves import Oval, Curve
import numpy as np
import matplotlib.pyplot as plt

from geometry import Projection, Segment, Point
from wurf import WURF


class Helper:

    def __init__(self):
        super(Helper, self).__init__()

    @staticmethod
    def solve_system(i, points):
        window = Helper.window(i, points)
        count = 5
        if len(window) != count:
            raise Exception("len(points) must be {}".format(count))
        A = np.asarray([
            np.asarray([p.x**2, p.y**2, p.x*p.y, p.x, p.y])
            for p in window
        ])
        B = np.asarray([1 for i in range(count)])
        return np.linalg.solve(A, B)

    @staticmethod
    def window(i, points, size=2):
        if i - size < 0:
            result = points[i - size:] + points[:i + size + 1]
        elif i + size >= len(points):
            result = points[i - size:] + points[:i + size + 1 - len(points)]
        else:
            result = points[i - size : i + size + 1]
        return result

    @staticmethod
    def calculate(points):
        result = []
        for i in range(0, len(points)):

            solution = Helper.solve_system(i, points)
            result.append(solution)
        return result

    @staticmethod
    def split(oval):
        up_points = oval.points[0:int(len(oval.points) / 2)]
        down_points = oval.points[int(len(oval.points) / 2) + 1:int(len(oval.points))]
        return up_points, down_points

    @staticmethod
    def draw_segment(s):
        x = [s.p1.x, s.p2.x]
        y = [s.p1.y, s.p2.y]
        plt.plot(x, y, [2])

    @staticmethod
    def draw_derivatives(points, step=10):
        for i in range(0, len(points) + 1, step):
            tg = Curve.tangent(*Helper.window(i, points))
            if tg:
                Helper.draw_segment(tg)

    @staticmethod
    def draw_sys_solutions(points, values):
        x = np.asarray([item.x for item in points])
        y = np.asarray([item[0] for item in values])
        plt.scatter(x, y, [2 for i in x])

    @staticmethod
    def find_conjugation_points(points, values, delta=0.1**5):
        result = []
        evaluate = lambda a, b: abs(a - b) < delta
        c1 = [item[0] for item in values]
        for i in range(len(points)):
            window = Helper.window(i, c1, 4)
            if (
                evaluate(window[0], window[1]) and
                evaluate(window[-2], window[-1]) and
                not evaluate(window[2], window[0]) and
                not evaluate(window[3], window[0]) and
                not evaluate(window[4], window[0]) #and
                # not evaluate(window[-4], window[-1])
            ):
                result.append(points[i - 1])
        return result

    @staticmethod
    def points_to_x_y(points):
        x = [item.x for item in points]
        y = [item.y for item in points]
        return x, y

    @staticmethod
    def get_inner_curve(oval, cross_point, wurf, step=1):
        curve = []
        for i in range(0, int(len(oval.points)), step):
            index = i + int(len(oval.points)/2)
            if index >= len(oval.points):
                index = i - int(len(oval.points)/2)
            cross = [oval.points[i], oval.points[index]]

            w = WURF.last_point(cross[0], cross_point, cross[1], wurf)
            curve.append(w.p3)
        return curve

    @staticmethod
    def wurf_mapping(first, second, main):
        result = []
        for i in range(0, len(first.points), 1):
            window = Helper.window(i, first.points, 2)
            tg = first.tangent(*window)
            [p2, p4] = second.cross_segment(tg)
            [p1, p5] = main.cross_segment(tg)
            p3 = first.points[i]
            [p1, p2, p3, p4, p5] = sorted([p1, p2, p3, p4, p5], key=lambda item: (item.x, item.y))
            w1 = WURF(p1, p2, p3, p4)
            w2 = WURF(p2, p3, p4, p5)
            result.append(Point(w1.value, w2.value))
            # x, y = Helper.points_to_x_y([p1, p2, p3, p4, p5])
            # plt.scatter(x, y)
        return result


    @staticmethod
    def main(oval, subplot=None):
        if subplot:
            plt.subplot(subplot)
        oval.draw()

        points = Helper.find_conjugation_points(oval.points, Helper.calculate(oval.points))
        cross_point = Segment.cross(Segment(points[0], points[2]), Segment(points[1], points[3]))
        x, y = Helper.points_to_x_y(points + [cross_point])
        plt.scatter(x, y)

        curve = Helper.get_inner_curve(oval, cross_point, 2, 1)
        first = Curve.from_points(tuple(curve))
        x, y = Helper.points_to_x_y(curve)
        plt.scatter(x, y, [2])

        curve = curve = Helper.get_inner_curve(oval, cross_point, 1.5, 1)
        second = Curve.from_points(tuple(curve))
        x, y = Helper.points_to_x_y(curve)
        plt.scatter(x, y, [2])

        wurf_map = Helper.wurf_mapping(first, second, oval)
        return wurf_map


if __name__ == '__main__':
    # w = WURF.last_point(Point(0, 0), Point(1, 0), Point(2, 0), 2)
    # print(w, w.value)

    oval = Oval(0.001)
    print(len(oval.points))
    wurf_map_1 = Helper.main(oval, 231)
    # Helper.main(oval, None)

    projected = Curve.from_points(Projection(
        1.5, 1, 0,
        1, 2, 0,
        0, 0.1,
    ).transform(oval.points))
    wurf_map_2 = Helper.main(projected, 232)

    projected_1 = Curve.from_points(Projection(
        1.5, 1, 0,
        1, 2, 0,
        0.2, 0.1,
    ).transform(oval.points))
    wurf_map_3 = Helper.main(projected_1, 233)

    plt.subplot(234)
    x, y = Helper.points_to_x_y(wurf_map_1)
    plt.scatter(x, y, [2])

    plt.subplot(235)
    x, y = Helper.points_to_x_y(wurf_map_2)
    plt.scatter(x, y, [2])

    plt.subplot(236)
    x, y = Helper.points_to_x_y(wurf_map_3)
    plt.scatter(x, y, [2])

    # Helper.draw_derivatives(oval.points, 10)
    # plt.subplot(222)
    # up, down = Helper.split(oval)
    # Helper.draw_sys_solutions(up, Helper.calculate(up))
    # plt.subplot(223)
    # Helper.draw_sys_solutions(down, Helper.calculate(down))

    plt.show()

