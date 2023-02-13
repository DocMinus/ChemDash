# Environment Setup

## Installation with Anaconda/Miniconda

Run the two commands from the root directory.

```shell
conda env create -f ./environment/conda.yaml
conda activate py311chemdash
```

## Installation with Pip
**Note:** Python >= 3.9 is required to use the [modern style](https://peps.python.org/pep-0585/) of type annotations.<br>
I prefer using 3.11 (due to increased performance over versions <=3.10)

Run the command from the root directory

```shell
python -m pip install -r ./environment/requirements.txt
```

## Problems with Windows install
dash-bio might give an install error due to parmed error. The latter requires MS Visual C++14.0 or greater.<br>
https://visualstudio.microsoft.com/visual-cpp-build-tools <br>
After that, rerun the above install again.


