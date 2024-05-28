from setuptools import setup, find_packages

setup(
    name='pydeseq',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'scipy',
        'matplotlib',
    ],
    entry_points={
        'console_scripts': [
            'pydeseq=pydeseq.pydeseq:main',
        ],
    },
)
