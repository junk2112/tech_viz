from geometry import Segment, Point
import numpy as np


class WURF:

    def __init__(self, p1, p2, p3, p4):
        self.p1, self.p2, self.p3, self.p4 = p1, p2, p3, p4
        self.a, self.b, self.c = Segment(p1, p2), Segment(p2, p3), Segment(p3, p4)

    def __str__(self):
        return '{}, {}, {}, {}'.format(self.p1, self.p2, self.p3, self.p4)

    def __repr__(self):
        return self.__str__()

    @property
    def value(self):
        la, lb, lc = self.a.len, self.b.len, self.c.len
        return (la + lb) * (lb + lc) / (lb * (la + lb + lc))

    @staticmethod
    def last_point_1(p1, p2, p4, wurf_value, delta=0.1**3):
        main = Segment(p1, p4)
        y = lambda x: main.k * x + main.b
        if p2.x > p4.x:
            delta *= -1
        current_x = p2.x
        wurfs = []
        while np.sign(delta)*current_x < p4.x:
            current_x += delta
            w = WURF(p1, p2, Point(current_x, y(current_x)), p4)
            # print(w.value)
            wurfs.append(w)
        tmp = [abs(item.value - wurf_value) for item in wurfs]
        return wurfs[tmp.index(min(tmp))]

    @staticmethod
    def last_point(p1, p2, p4, wurf_value):
        main = Segment(p1, p4)
        y = lambda x: main.k * x + main.b
        a = Segment(p1, p2)
        b_c = Segment(p2, p4)
        l = a.len / (((wurf_value * (a.len + b_c.len)) / b_c.len) - 1)
        l = l**2
        ua = 1 + (main.k**2)
        ub = 2 * (main.k * main.b - p2.x - p2.y * main.k)
        uc = (p2.x**2) + (p2.y**2) + (main.b**2) - 2*p2.y*main.b - l
        D = (ub ** 2) - (4 * ua * uc)
        if D < 0:
            return None
        x1 = (-ub + D ** 0.5) / (2 * ua)
        x2 = (-ub - D ** 0.5) / (2 * ua)
        w1, w2 = WURF(p1, p2, Point(x1, y(x1)), p4), WURF(p1, p2, Point(x2, y(x2)), p4)
        result = w1
        if (
            (abs(w2.value - wurf_value) < abs(w1.value - wurf_value))
        ):
            result = w2
        return result
