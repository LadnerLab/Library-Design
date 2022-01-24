---
layout: default
title: Installation
permalink: /installation/
---

## Installation

Installation is based off of each independent module
  - C-KmerOligo
  - ...

### Installation of C-KmerOligo

#### Dependencies

- Linux/OSX Operating System, all other files/libraries are included
- GCC

#### Installation
##### MacOS
You must have make and gcc installed in order to build this program, this
cab be achieved with
```
xcode-select --install
```
from the command line or through Xcode itself. Once you have these installed,
follow the instructions for Linux

##### Linux
```
git clone https://github.com/LadnerLab/c-KmerOligo.git
cd C-KmerOligo
make optimized
```
