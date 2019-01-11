#!/bin/bash

# Runs allostery ID
# Mounts a directory in the container that can be viewed in the data file of the working directory
sudo docker run -it -v $(pwd)/data:/home/luser/data allostery-id