from setuptools import setup

setup(
    name="samUtilities",
    version="0.9.0",
    packages=["samUtilities","samUtilities.betterHits"],
    test_suite="tests",
    url="https://github.com/gdbzork/samUtilities",
    license="MIT ",
    author="Gord Brown",
    author_email="gdbzork@gmail.com",
    description="Miscellaneous SAM/BAM Utilities",
    scripts=["bin/reportBetterHits"]
)
