project('purelib-and-platlib', 'c')

py = import('python').find_installation(pure: false)

py.extension_module(
    'example_module',
    'example_module.c',
    install: true,
)
