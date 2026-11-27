FROM python:3.13-slim
WORKDIR /opt/spatial-audit
COPY . .
RUN pip install --no-cache-dir ".[agent,dev,archives]"
ENTRYPOINT ["spatial-audit"]
CMD ["demo", "--out", "/tmp/spatial-audit-results"]
