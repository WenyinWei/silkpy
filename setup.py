import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="silkpy",
    version="0.0.1",
    author="Wenyin Wei",
    author_email="wenyin.wei@ipp.ac.cn",
    description="Silkpy is planned to be a math package for curves and surfaces which are classical objects of differential geometry. Silk, in Chinese, is a kind of soft clothe material with intervined perpendicular lines of longitude and latitude. The lines are named '丝 si' while the clothe is named '绸 chou', which are metaphors for the 'curve' and 'surface' in differential geometry.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/WenyinWei/silkpy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)