from setuptools import find_packages, setup

setup(
    name="gem-utilities",  # Changed from 'GEM-utilities' to match package folder name
    version="0.1dev",
    packages=find_packages(),  # Automatically finds all sub-packages
    license="MIT",
    description="Utilities for GEM models",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    install_requires=[],  # Add required dependencies here
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
