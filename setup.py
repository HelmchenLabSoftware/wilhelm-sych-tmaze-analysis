import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="chernysheva-tmaze-analysis", # Replace with your own username
    version="0.0.1",
    author="Aleksejs Fomins",
    author_email="fomins.aleksejs@gmail.com",
    description="Code for PFC memory maintainance publication by Sych et al 2023",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aleksejs-fomins/chernysheva-tmaze-analysis",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)
