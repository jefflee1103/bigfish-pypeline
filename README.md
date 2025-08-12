# Big-FISH FISH Image Quantification and Analysis Pypeline

Pre-process and quantify single-molecule Fluorescence in situ Hybridisation images using a python package Big-FISH: https://github.com/fish-quant/big-fish.


## Installation

### Install Pixi

Pixi is a package/environment manager. Install using the official guide: https://pixi.sh/dev/installation/ .

#### Linux & MacOS
    
    curl -fsSL https://pixi.sh/install.sh | sh

#### Windows

    powershell -ExecutionPolicy ByPass -c "irm -useb https://pixi.sh/install.ps1 | iex"


### Setup pixi environment for the BigFISH pypeline

Clone or download the contents of this repo and change directory to the root directory of this repo.

    cd /path/to/repo/root/directory

    pixi install


### Running the pyeline

    # (option 1) pixi default
    pixi run <your task>

    # (option 2) conda style
    pixi shell
    <your task>
    exit # for exiting the envornment






