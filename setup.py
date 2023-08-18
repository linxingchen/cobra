from setuptools import setup, find_packages

setup(
    name="cobra",
    version="1.1.54",
    description="a bioinformatics tool to get higher quality viral genomes assembled from metagenomes of short paired-end reads",
    long_description=open("README.md", "r").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/linxingchen/cobra.github.io",
    author="Linxing Chen",
    author_email="linkingchan@gmail.com",
    license="MIT",  # or other licenses you prefer
    classifiers=[
        "Development Status :: 3 - Alpha",  # or "5 - Production/Stable" or another status
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: MIT License",  # matching the license from above
        "Programming Language :: Python :: 3",  # or "2.7", "3.6", etc.
    ],
    python_requires=">=3.7",  # your package's Python version requirements
)
