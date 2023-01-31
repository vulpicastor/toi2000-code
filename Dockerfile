FROM continuumio/miniconda3:4.12.0

RUN conda update -y --update-all
RUN conda install astropy matplotlib pymc3 pytables
RUN conda install -c conda-forge celerite2 corner cxx-compiler 'exoplanet-core<0.2' exoplanet jupyterlab pydl pymc3-ext
RUN mkdir -p /opt/git/toi2000-code

EXPOSE 8888

ENTRYPOINT jupyter-lab --notebook-dir=/opt/git/toi2000-code --ip='*' --port=8888 --no-browser --allow-root
