import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='measure_process',
    version='0.0.1',
    description='Testing installation of Package',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/opietx/Measure_process',
    license='MIT',
    packages=['measure_process'],
)