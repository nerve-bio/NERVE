FROM francecosta/nerve:v0.0.1 AS final
RUN pip install epitopepredict
RUN pip install pandas==2.0.2

COPY ./ /usr/nerve_python/NERVE
WORKDIR /workdir
EXPOSE 8880
RUN chmod +x /usr/nerve_python/NERVE/code/NERVE.py

ENTRYPOINT ["/usr/nerve_python/NERVE/code/NERVE.py"]
