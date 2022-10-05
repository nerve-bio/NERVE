### Run the image as a container with a local folder as "work" folder,
### and one port exposed for access to Jupyter notebook (with browser)
docker run -d -p 8880:8880 -v $HOME/code/nerve:/my_data nerve
