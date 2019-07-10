FROM pinellolab/stream

RUN conda install -c anaconda pyyaml -y

COPY main.py runner.sh definition.yml /code/

RUN chmod +x /code/*

ENTRYPOINT ["/code/runner.sh"]
