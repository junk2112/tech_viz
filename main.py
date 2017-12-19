import sys
sys.setrecursionlimit(10000)
import itertools
from functools import reduce
from geometry import Frechet

from curves import Oval, Curve, Ellipse, Circle
import numpy as np
import matplotlib.pyplot as plt
from pyflann import *

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
        # print(A)
        # print(B)
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
            try:
                solution = Helper.solve_system(i, points)
                result.append(solution)
            except:
                sys.exit()
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
    def find_conjugation_points(points, values, delta=0.1**4):
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
        if len(result) < 4:
            raise Exception('Found only {} conjugation points'.format(len(result)))
        if len(result) > 4:
            distances = sorted(
                [
                    (
                        comb,
                        # reduce(lambda a, b: a * b, [Segment(comb[0], comb[i]).len for i in range(1, len(comb))]) *
                        reduce(lambda a, b: a * b, [abs(points.index(comb[i]) - points.index(comb[i-1])) for i in range(len(comb))])
                    )
                    for comb in itertools.combinations(result, 4)
                ],
                reverse=True,
                key=lambda item: item[1]
            )
            result = list(distances[0][0])
            result = sorted(
                [
                    (item, points.index(item))
                    for item in result
                ],
                key=lambda item: item[1]
            )
            result = [item[0] for item in result]
        return result

    @staticmethod
    def points_to_x_y(points):
        x = [item.x for item in points]
        y = [item.y for item in points]
        return x, y

    @staticmethod
    def get_inner_curve(oval, cross_point, wurf, step=1):
        def get_p3(p1, p2):
            cross = oval.cross_segment(Segment(p1, p2))
            w = WURF.last_point(cross[0], cross_point, cross[1], wurf)
            if not w:
                # print('getting wurf: D < 0')
                return None
                raise Exception('getting wurf: D < 0')
            return w.p3

        result = []
        l_2 = int(len(oval.points)/2)
        for i in range(0, len(oval.points), step):
            index = i + l_2
            if index >= len(oval.points):
                index = i - l_2
            p3 = get_p3(oval.points[index], oval.points[i])
            result.append(p3)
        return result

    @staticmethod
    def wurf_mapping(first, second, main):
        result = []
        for i in range(0, len(first.points), 1):
            window = Helper.window(i, first.points, 2)
            tg = first.tangent(*window)
            try:
                [p2, p4] = second.cross_segment(tg)
                [p1, p5] = main.cross_segment(tg)
            except:
                continue
            else:
                p3 = first.points[i]
                [p1, p2, p3, p4, p5] = sorted(
                    [p1, p2, p3, p4, p5],
                    key=lambda item: (item.x*tg.k, item.y)
                )
                w1 = WURF(p1, p2, p3, p5)
                w2 = WURF(p1, p2, p4, p5)
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
        print(cross_point)
        x, y = Helper.points_to_x_y(points + [cross_point])
        # x, y = Helper.points_to_x_y(points)

        plt.scatter(x, y)

        curve = Helper.get_inner_curve(oval, cross_point, 2, 1)
        first = Curve.from_points(tuple(curve), step)
        x, y = Helper.points_to_x_y(curve)
        plt.scatter(x, y, [2])

        curve = curve = Helper.get_inner_curve(oval, cross_point, 1.5, 1)
        second = Curve.from_points(tuple(curve), step)
        x, y = Helper.points_to_x_y(curve)
        plt.scatter(x, y, [2])

        wurf_map = Helper.wurf_mapping(first, second, oval)
        result = []
        for i in range(len(wurf_map)):
            nearest, dist = FLANN().nn(
                np.asarray([(item.x, item.y) for item in wurf_map if item != wurf_map[i]]),
                np.asarray([(item.x, item.y) for item in wurf_map if item == wurf_map[i]]),
                1, algorithm="kmeans",
            )
            # print(dist)
            if dist[0] < 0.01:
                result.append(wurf_map[i])
            # else:
            #     print('here')
        return result

    @staticmethod
    def main_main(oval, number):
        number *= 4
        # print(len(oval.points))
        # wurf_map_1 = Helper.main(oval, 241 + number)
        plt.subplot(241 + number)
        oval.draw()

        projected_1 = Curve.from_proj(Projection(
            1.5, 1, 0,
            1, 2, 0,
            0, 0.2,
        ), oval.points[0::2], step)
        print(len(projected_1.points))
        wurf_map_2 = Helper.main(projected_1, 242 + number)
        #
        projected_2 = Curve.from_proj(Projection(
            1.5, 1, 0,
            1, 2, 0,
            0.2, 0.1,
        ), oval.points[1::2], step)
        print(len(projected_2.points))
        wurf_map_3 = Helper.main(projected_2, 243 + number)

        # projected_2 = Curve.from_proj(Projection(
        #     1.5, 1, 1,
        #     1, 2, 1,
        #     0.2, 0.1,
        # ), oval.points[1::2], step)
        # wurf_map_3 = Helper.main(projected_2, 243 + number)

        # plt.subplot(244)
        # x, y = Helper.points_to_x_y(wurf_map_1)
        # plt.scatter(x, y, [2])
        #
        plt.subplot(244 + number)
        x, y = Helper.points_to_x_y(wurf_map_2)
        plt.scatter(x, y, [2])
        #
        plt.subplot(244 + number)
        x, y = Helper.points_to_x_y(wurf_map_3)
        plt.scatter(x, y, [2])

        print('Frechet distance is')
        print(Frechet.dist(wurf_map_2, wurf_map_3))
        return [wurf_map_2, wurf_map_3]

step = 0.0001

if __name__ == '__main__':


    e1 = Ellipse(5, 1, step)
    e2 = Circle(1, 0, 0, step)
    center = Point(2, 0)

    oval = Oval(e1, e2, center, step)
    m1, m2 = Helper.main_main(oval, 0)

    e1 = Ellipse(6, 1, step)
    e2 = Ellipse(2, 1, step)
    center = Point(2, 0)

    oval = Oval(e1, e2, center, step)
    m3, m4 = Helper.main_main(oval, 1)

    # print(Frechet.dist(m1, m3))
    # print(Frechet.dist(m2, m4))

    plt.show()

