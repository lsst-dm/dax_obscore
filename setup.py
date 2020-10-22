from setuptools import setup

setup(
    name='obscore_generator',
    package_dir={'lsst': 'python/lsst'},
    package_data={'lsst': ['dax/obscore_generator/config/*']},
    packages=['lsst', 'lsst.dax', 'lsst.dax.obscore_generator'],
    zip_safe=False,
    use_scm_version={'version_scheme': 'post-release'},
    setup_requires=['setuptools_scm'],
)
