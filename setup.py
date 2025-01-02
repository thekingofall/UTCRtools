
from setuptools import setup, find_packages

setup(
    name='UTCRtools',
    version='0.2.0',
    description='A TCR data processing toolkit',
    author='Your Name',
    author_email='your_email@example.com',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'biopython',
        # 如果脚本还需要别的库，请在这里添加，如 'numpy', 'configparser' 等
    ],
    entry_points={
        'console_scripts': [
            # 这个 utcrtools 是安装后命令行调用的名字
            # 指向 utcrtools/main.py 里的 main() 函数
            'utcrtools=utcrtools.main:main',
        ],
    },
    python_requires='>=3.6',
)
