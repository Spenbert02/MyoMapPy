
py.install_sources(
    [
        '__init__.py',
    ],
    subdir: 'myomappy',
)

subdir('c_utils')  # this needs to be first, so tractography can include it
subdir('tractography')
