"""
4th order Runge Kutta Tracking
"""

from myomappy.rungekutta import _runge_kutta_4_calcs

def rk4Track():
    print("4th order runge kutta track")
    print(_runge_kutta_4_calcs.foo())
