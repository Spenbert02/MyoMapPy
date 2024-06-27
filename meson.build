project('purelib-and-platlib', 'c')

py = import('python').find_installation(pure: false)

py.extension_module(
    'runge_kutta_4_calcs',
    'myomappy/rungekutta/runge_kutta_4_calcs.c',
    install: true,
)
