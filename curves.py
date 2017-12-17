from geometry import Point, Segment
import numpy as np
import matplotlib.pyplot as plt
import math
import sys


class Curve:

    def __init__(self, step):
        super(Curve, self).__init__()
        self.points = []
        self.step = step
        self.default_range = (-10, 10)

    def add_points(self, p1, p2):
        if (
                Segment(p1, p2).len > self.step
                and p1.parent is not None
                and p1.parent == p2.parent
        ):
            x = (p1.x + p2.x)/2
            y = p1.parent.y(x)
            y = y[0] if abs(y[0] - p1.y) < abs(y[1] - p1.y) else y[1]
            point = Point(x, y, self)
            for p in self.add_points(p1, point):
                yield p
            yield point
            for p in self.add_points(point, p2):
                yield p

    def handle_point(self, p1, p2):
        if p2.x is None or p2.y is None:
            return []
        new = list(self.add_points(p1, p2))
        return new + [p2]

    def get_points(self):
        a1, a2 = self.default_range
        r1, r2 = [], []
        current = a1
        while current <= a2:
            y = self.y(current)
            p1 = Point(current, y[0], self)
            p2 = Point(current, y[1], self)
            # if p1.y is not None:
            #     if len(r1) == 0:
            #         r1.append(p1)
            #     else:
            #         r1 += self.handle_point(r1[-1], p1)
            # if p2.y is not None:
            #     if len(r2) == 0:
            #         r2.append(p2)
            #     else:
            #         r2 += self.handle_point(r2[-1], p2)
            r1.append(p1)
            r2.append(p2)
            current += self.step
        f_l = lambda p: (
                (p.x is not None) and
                (p.y is not None)
            )
        return list(filter(f_l, r1)) + list(filter(f_l, r2))

    def draw(self):
        x = np.asarray([p.x for p in self.points])
        y = np.asarray([p.y for p in self.points])
        plt.scatter(x, y, [2 for i in x])
        # plt.show()

    @staticmethod
    def derivative(a1, a2, c1, c2):
        return (
            (a1.y + a2.y - c1.y - c2.y) /
            (a1.x + a2.x - c1.x - c2.x)
        )

    @staticmethod
    def tangent(a1, a2, current, c1, c2):
        try:
            k = Curve.derivative(a1, a2, c1, c2)
        except ZeroDivisionError:
            return None
        else:
            b = current.y - k * current.x
            return Segment.from_line_and_point(k, b, current)

    def cross_segment(self, segment):
        result = []
        y = lambda x: segment.k * x + segment.b
        cache = {"{}-{}".format(p.x, p.y): Segment(p, Point(p.x, y(p.x))).len for p in self.points}
        def min_p(points, prev=None, delta=1):
            r = []
            if not prev:
                dists = [cache["{}-{}".format(p.x, p.y)] for p in points]
            else:
                dists = [cache["{}-{}".format(p.x, p.y)] / Segment(p, prev).len for p in points]
            m_v = min(dists)
            if m_v < delta:
                r.append(points[dists.index(m_v)])
            return r
        result += min_p(self.points)
        if result:
            index = self.points.index(result[0])
            points = self.points[:index] + self.points[index+1:]
            result += min_p(points, prev=result[-1])
        return list(sorted(result, key=lambda item: Segment(item, segment.p1).len))

    def y(self, x):
        raise NotImplementedError()

    @staticmethod
    def from_points(points, step):
        r = Curve(step)
        r.points = points
        return r

    @staticmethod
    def from_proj(projection, parent_points, step):
        points = list(projection.transform(parent_points))
        result = Curve(step)

        # resampled = [points[0]]
        # for i in range(1, len(points)):
        #     if Segment(resampled[-1], points[i]).len > step:
        #         resampled += list(projection.transform(
        #             result.handle_point(parent_points[i-1], parent_points[i])
        #         ))
        # result.points = resampled
        # result.points = result.filter_nearest(result.points, result.step)

        result.points = points
        return result

    def filter_nearest(self, points, delta):
        l_2 = int(len(points)/2)
        r1, r2 = [points[0]], [points[-1]]
        for i in range(1, l_2):
            c1, c2 = points[i], points[len(points) - i]
            l1, l2 = r1[-1], r2[-1]
            if Segment(c1, l1).len**2 > delta and Segment(c2, l2).len**2 > delta:
                r1.append(c1)
                r2.append(c2)
        return r1 + r2[::-1]

class Ellipse(Curve):

    def __init__(self, a, b, step, offset_x=0, offset_y=0):
        super(Ellipse, self).__init__(step)
        self.a, self.b = a, b
        self.offset_x, self.offset_y = offset_x, offset_y
        self.a_2, self.b_2 = self.a ** 2, self.b ** 2
        self.b_der_a = b / a
        self.excent = (self.b_2 / self.a_2) ** 0.5
        self.foc_par = self.b_2 / self.a
        self.points = self.get_points()

    def y(self, x):
        pre_result = self.pre_y(x)
        if isinstance(pre_result, float):
            return pre_result + self.offset_y, -pre_result + self.offset_y
        else:
            return None, None

    def pre_y(self, x):
        return self.b_der_a * (self.a_2 - (x - self.offset_x) ** 2) ** 0.5

    def radius(self, angle):
        return self.foc_par/(1 + self.excent * math.cos(angle))

class Circle(Ellipse):

    def __init__(self, r, offset_x, offset_y, step):
        super(Circle, self).__init__(r, r, step, offset_x, offset_y)


class Oval(Curve):

    def __init__(self, e1, e2, center, step):
        super(Oval, self).__init__(step)
        points = []
        points += list(filter(
            lambda p: (
                p.x < e2.offset_x and
                p.y > e2.offset_y
            ), e2.get_points()
        ))
        points += list(filter(
            lambda p: (
                p.y >= e1.offset_y and
                p.x >= e1.offset_x
            ), e1.get_points()
        ))
        points += [p.sim(center) for p in points]
        y_pos = list(filter(lambda p: p.y >= 0, points))
        y_neg = list(filter(lambda p: p.y < 0, points))
        y_pos = sorted(y_pos, key=lambda p: p.p_angle)
        y_neg = sorted(y_neg, key=lambda p: p.p_angle, reverse=True)
        y_pos, y_neg = self.filter_nearest(y_pos, self.step), self.filter_nearest(y_neg, self.step)
        self.points = y_pos + y_neg
