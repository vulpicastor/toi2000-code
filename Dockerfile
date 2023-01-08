FROM continuumio/miniconda3:4.12.0

RUN conda update -y --update-all
RUN conda install pymc3
RUN conda install matplotlib
RUN conda install astropy
RUN conda install -c conda-forge celerite2
RUN conda install -c conda-forge corner
RUN conda install -c conda-forge cxx-compiler
RUN conda install -c conda-forge 'exoplanet-core<0.2'
RUN conda install -c conda-forge exoplanet
RUN conda install -c conda-forge jupyterlab
RUN conda install -c conda-forge pymc3-ext
RUN mkdir -p /opt/git/toi2000-code

EXPOSE 8888

ENTRYPOINT jupyter-lab --notebook-dir=/opt/git/toi2000-code --ip='*' --port=8888 --no-browser --allow-root
