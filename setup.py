from setuptools import setup, find_packages

setup(
    name="ttpsolver",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        'numpy==1.22.1',
        'matplotlib==3.5.1',
        'pandas==1.4.1',
        'astropy==5.0.1',
        'astroplan==0.8',
        'plotly==5.17.0',
        'jupyter_client==7.1.2',
        'ipykernel==6.9.1',
        'imageio==2.24.0',
        'kaleido==0.2.1',
        'gurobipy'
    ],
    description="The Traveling Telescope Problem Solver, orders a set of targets in the optimal slew path.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Luke Handley",
    author_email="lhandley@caltech.edu",
    url="https://github.com/lukehandley/ttp/",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)