FROM francecosta/nerve:v0.0.1 AS final
RUN pip install epitopepredict

WORKDIR /workdir
RUN pip install epitopepredict
EXPOSE 8880
RUN chmod +x /usr/nerve_python/NERVE/code/NERVE.py

#CMD jupyter notebook --notebook-dir=/workdir --ip='0.0.0.0' --port='8880' --allow-root --NotebookApp.token='' --NotebookApp.password=''
ENTRYPOINT ["/usr/nerve_python/NERVE/code/NERVE.py"]
