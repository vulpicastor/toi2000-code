# Build instructions

If the user has an existing Anaconda installation, it suffices to create a new
environment from the supplied `environment.yml` to reproduce the environment
used in the preparation of the published paper.

If the user wishes to use Docker or Podman to create a self-contained
installation for running the code in this repository, brief instructions are
provided below. This is considered an expert option.

Podman (on Fedora 37) was used to run the container by the author, but in
theory the following instructions also apply to Docker by substituting the
`podman` command for `docker`.

## Creating and running the production container

First, build the production container using the following command with the
present working directory set to the root of this repository.
```
podman build -f Dockerfile
```
The command should output a hash of the built container image. It can then be
run using the following command, substituting `$IMAGEHASH` with the actual
hash of the image.
```
podman run -i -t -d --name toi2000 -p 8888:8888 --volume "$PWD":/opt/git/toi2000-code $IMAGEHASH
```
This recreates a new container named `toi2000` and binds the exposed port 8888
to localhost. The URL for accessing the Jupyter Lab instance can then be found
by the following command.
```
podman exec toi2000 jupyter-server list
```
Direct your browser to the URL listed to access the Jupyter Lab instance
running inside of the container.


## Creating and running the development container

The `environment.yml` file was generated from a development version of the
container. If you wish to recreate the development container with the latest
version of packages, you can substitute the build command with the following.
```
podman build -f Dockerfile.dev
```
The instructions for running the development container is identical to those
for the production container above.


## Special consideration for SELinux

On distributions that use SELinux for app isolation, like the desktop editions
of Fedora, it is necessary to relabel the SELinux constext of the files in
this repository in order to make them writable from within the container. This
can be accomplished by the following command.
```
chcon -R -t container_file_t .
```
