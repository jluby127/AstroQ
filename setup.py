from setuptools import setup, find_packages

setup(
    name="kpfcc",
    version="0.1",
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
    python_requires="==3.9",
)
