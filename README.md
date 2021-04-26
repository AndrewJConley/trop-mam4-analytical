# trop-mam4-analytical
Analytical solution of trop-mam4 gas-phase solver

# Build and run in a container
You can run the tests and generate plots of the results that you can view locally by first creating a folder on your local machine that you can use to exchange data with the container:

```
mkdir /path/to/my/local/folder
```

Then, build the image and run it in a container, mounting your local folder in the container (note that absolute paths must be used for the `-v` arguments):

```
docker build -t trop-mam4-test .
docker run -v /path/to/my/local/folder:/local_folder -it trop-mam4-test bash
```

Inside the container, run the tests, plot the results, and copy the plots to the shared folder:

```
cd /build
make test
cd test
gnuplot plot.conf
cp out/* /local_folder/
```

You can then navigate to the folder on your local machine to view the plots
