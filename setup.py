from setuptools import setup, find_packages

setup(
    name='colbuilder',
    version='0.1.0',
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "click",
        "pydantic>=2.0",
        "pyyaml",
        "colorama"
    ],
    entry_points={
        'console_scripts': [
            'colbuilder=colbuilder.colbuilder:main',
        ],
    },
)
