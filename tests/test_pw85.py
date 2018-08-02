import pytest

import pw85

def test_foo():
    q = pw85.spheroid(1.0, 0.1, None)
    for qi in q:
        print(qi)
