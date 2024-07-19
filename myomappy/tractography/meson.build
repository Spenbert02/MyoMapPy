# shared_library(  # cant get around segmentation faults. this library is nonfunctional as of 7/19/2024
#     '_tractography',
#     'tractography.c',
#     include_directories: [
#         'tractography_sources',
#         c_utils
#     ],
#     install: true,
#     install_dir: py.get_install_dir() / 'myomappy/tractography'
# )

py.install_sources(
    [
        '__init__.py',
        '_utils.py',
        'rungekutta.py',
        'euler.py'
    ],
    subdir: 'myomappy/tractography',
)
