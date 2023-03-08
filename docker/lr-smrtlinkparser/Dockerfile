FROM python:3.9.16-slim-buster

RUN apt-get update \
  && apt-get install -y build-essential python3-lxml --no-install-recommends \
  && pip install xmltodict==0.12.0 \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/doc && rm -rf /usr/share/man \
  && apt-get purge -y --auto-remove build-essential \
  && apt-get clean

ENV PYTHONUNBUFFERED="true"
