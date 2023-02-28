FROM continuumio/miniconda3:4.12.0

RUN conda env create -f environment.yml
RUN mkdir -p /opt/git/toi2000-code

EXPOSE 8888

ENTRYPOINT jupyter-lab --notebook-dir=/opt/git/toi2000-code --ip='*' --port=8888 --no-browser --allow-root
