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

    def get_points(self):
        a1, a2 = self.default_range
        result = []
        current = a1
        while current <= a2:
            y = self.y(current)
            result.append(Point(current, y[0]))
            result.append(Point(current, y[1]))
            current += self.step
        result = list(filter(
            lambda p: (
                (p.x is not None) and
                (p.y is not None)
            ),
            result
        ))
        return result

    def draw(self):
        x = np.asarray([p.x for p in self.points])
        y = np.asarray([p.y for p in self.points])
        plt.scatter(x, y, [2 for i in x])
        # plt.show()

    @staticmethod
    def derivative(a1, a2, c1, c2):
        return (a1.y + a2.y - c1.y - c2.y) / (a1.x + a2.x - c1.x - c2.x)

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
        def min_p(points, prev=None, delta=1):
            r = []
            if not prev:
                dists = [Segment(p, Point(p.x, y(p.x))).len for p in points]
            else:
                dists = [Segment(p, Point(p.x, y(p.x))).len / Segment(p, prev).len for p in points]
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
        # return result

    def y(self, x):
        raise NotImplementedError()

    @staticmethod
    def from_points(points):
        result = Curve(1)
        result.points = list(points)
        return result

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


class Pear(Ellipse):

    def __init__(self, a, b, step, pow_x=1., offset_x=0, offset_y=0):
        self.pow_x = pow_x
        super(Pear, self).__init__(a, b, step, offset_x, offset_y)

    def pre_y(self, x):
        return self.b_der_a * (self.a_2 - (x**self.pow_x - self.offset_x) ** 2) ** 0.5

class Circle(Ellipse):

    def __init__(self, r, offset_x, offset_y, step):
        super(Circle, self).__init__(r, r, step, offset_x, offset_y)


class Oval(Curve):

    def __init__(self, step):
        super(Oval, self).__init__(step)
        ellipse = Ellipse(5, 1, self.step)
        circle = Circle(1, 0, 0, self.step)
        points = []
        points += list(filter(
            lambda p: (
                p.x < 0 and p.y > 0
            ), circle.get_points()
        ))
        points += list(filter(
            lambda p: (
                p.y >= ellipse.offset_y and
                p.x >= ellipse.offset_x
            ), ellipse.get_points()
        ))
        points += [p.sim(Point(2, 0)) for p in points]
        y_pos = list(filter(lambda p: p.y >= 0, points))
        y_neg = list(filter(lambda p: p.y < 0, points))
        y_pos = sorted(y_pos, key=lambda p: p.p_angle)
        y_neg = sorted(y_neg, key=lambda p: p.p_angle, reverse=True)
        y_pos, y_neg = self.filter_nearest(y_pos, self.step), self.filter_nearest(y_neg, self.step)
        self.points = y_pos + y_neg

    def filter_nearest(self, points, delta):
        result = [points[0]]
        for i in range(1, len(points)):
            current = points[i]
            last = result[-1]
            if (current.x - last.x)**2 + (current.y - last.y)**2 > delta:
                result.append(current)
        return result




if __name__ == '__main__':
    # e = Pear(1, 1, step=0.0001, pow_x=1.4, offset_x=1)
    # e.draw()
    Oval(0.001).draw()
    # Ellipse(1, 1, 0.01).draw()
    # Circle(1, 0, 0, 0.01).draw()