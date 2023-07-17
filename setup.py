from setuptools import setup

with open('./requirements.txt') as f:
    INSTALL_REQUIRES = f.read().splitlines()

setup(
    name='river_route',
    version='0.0.1',
    description='Perform river routing computations on large river networks',
    author='Riley Hales PhD',
    install_requires=INSTALL_REQUIRES
)
