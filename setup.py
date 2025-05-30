from setuptools import setup, find_packages
import re


def get_property(prop, project):
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop),
                       open(project + '/__init__.py').read())
    return result.group(1)

setup(
    name="kpfcc",
    version=get_property('__version__', 'kpfcc'),
    packages=find_packages(),
    install_requires=[],
    description="The Keck Planet Finder Community Cadence (KPFCC) auto-scheduler software.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Jack Lubin",
    author_email="jblubin@ucla.edu",
    url="https://github.com/jluby127/optimalAllocation/",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points = {
        "console_scripts": ['astroq=kpfcc.cli:main']
    },
    # python_requires="==3.9",
)
