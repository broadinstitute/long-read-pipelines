IMAGE_NAME = lr-splicedbam2gff
# Pinned to version of splicedbam2gff:
VERSION = 1.3

TAG1 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):latest

all: | build push

build:
    docker build -t $(TAG1) -t $(TAG2) .

build_no_cache:
    docker build --no-cache -t $(TAG1) -t $(TAG2) .

push:
    docker push $(TAG1)
    docker push $(TAG2)
