import math

class Point:

    def __init__(self, x, y):
        self.x, self.y = x, y

    def __str__(self):
        return "({}, {})".format(self.x, self.y)

    def __repr__(self):
        return self.__str__()

    def sim(self, center):
        return Point(2 * center.x - self.x, 2 * center.y - self.y)

    @property
    def p_radius(self):
        return (self.x**2 + self.y**2) ** 0.5

    @property
    def p_angle(self):
        return math.acos(self.x/self.p_radius)

    @staticmethod
    def from_polar(radius, angle):
        return Point(radius * math.cos(angle), radius * math.sin(angle))


class Segment:

    def __init__(self, p1, p2):
        self.p1, self.p2 = p1, p2
        try:
            self.k = (self.p1.y - self.p2.y) / (self.p1.x - self.p2.x)
            self.b = self.p1.y - self.k * self.p1.x
        except:
            self.k, self.b = None, None

    def __str__(self):
        return "[{}, {}]".format(self.p1, self.p2)

    @staticmethod
    def cross(s1, s2):
        def validate(point):
            return (
                min(s1.p1.x, s1.p2.x)
                <= point.x
                <= max(s1.p1.x, s1.p2.x)
            )

        A1 = s1.p1.y - s1.p2.y
        B1 = s1.p2.x - s1.p1.x
        C1 = s1.p1.x * s1.p2.y - s1.p2.x * s1.p1.y
        A2 = s2.p1.y - s2.p2.y
        B2 = s2.p2.x - s2.p1.x
        C2 = s2.p1.x * s2.p2.y - s2.p2.x * s2.p1.y

        result = None
        if B1 * A2 - B2 * A1 and A1:
            y = (C2 * A1 - C1 * A2) / (B1 * A2 - B2 * A1)
            x = (-C1 - B1 * y) / A1
            result = Point(x, y)
        elif B1 * A2 - B2 * A1 and A2:
            y = (C2 * A1 - C1 * A2) / (B1 * A2 - B2 * A1)
            x = (-C2 - B2 * y) / A2
            result = Point(x, y)
        # if result and not validate(result):
        #     result = None
        return result

    @staticmethod
    def from_line_and_point(k, b, point, delta=0.1):
        x1, x2 = point.x - delta, point.x + delta
        y = lambda x: k*x + b
        return Segment(Point(x1, y(x1)), Point(x2, y(x2)))

    @property
    def len(self):
        return (
            (self.p1.x - self.p2.x) ** 2 +
            (self.p1.y - self.p2.y) ** 2
        ) ** 0.5


class Projection:

    def __init__(self, a1, a2, a3, b1, b2, b3, c1, c2):
        self.a1, self.a2, self.a3 = a1, a2, a3
        self.b1, self.b2, self.b3 = b1, b2, b3
        self.c1, self.c2 = c1, c2

    def mapper(self, p):
        return Point(
            (self.a1 * p.x + self.a2 * p.y + self.a3) / (self.c1 * p.x + self.c2 * p.y + 1),
            (self.b1 * p.x + self.b2 * p.y + self.b3) / (self.c1 * p.x + self.c2 * p.y + 1),
        )

    def transform(self, points):
        result = []
        for p in points:
            result.append(self.mapper(p))
        return result
