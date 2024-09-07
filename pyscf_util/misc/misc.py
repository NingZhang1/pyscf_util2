def _combine2(a, b):
    if a > b:
        return a * (a + 1) // 2 + b
    else:
        return b * (b + 1) // 2 + a


def _combine4(a, b, c, d):
    return _combine2(_combine2(a, b), _combine2(c, d))
