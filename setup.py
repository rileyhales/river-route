import re

from setuptools import setup, find_packages

NAME = 'river_route'
DESCRIPTION = 'Perform river routing computations on large river networks'
URL = 'https://github.com/rileyhales/river-route'
AUTHOR = 'Riley Hales PhD'
REQUIRES_PYTHON = '>=3.12.0'

with open(f'./{NAME}/__metadata__.py') as f:
    version_pattern = r'__version__ = [\'"](\d+\.\d+\.\d+)[\'"]'
    VERSION = re.search(version_pattern, f.read()).group(1)
    print(f'Version: {VERSION}')

with open('./requirements.txt') as f:
    INSTALL_REQUIRES = f.read().splitlines()

with open('./docs/requirements.txt') as f:
    DOCS_REQUIRES = f.read().splitlines()

TEST_REQUIRES = ['pytest', ]
APP_REQUIRES = [f'river-route-app~={VERSION}', ]

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author=AUTHOR,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    install_requires=INSTALL_REQUIRES,
    extras_require={
        'all': TEST_REQUIRES + DOCS_REQUIRES + APP_REQUIRES,
        'test': TEST_REQUIRES,
        'docs': DOCS_REQUIRES,
        'app': APP_REQUIRES,
    },
    include_package_data=False,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: 3.14',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: GIS',
        'Topic :: Scientific/Engineering :: Hydrology',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
    entry_points={
        'console_scripts': ['rr=river_route._cli:main', ]
    },
    options={}
)
