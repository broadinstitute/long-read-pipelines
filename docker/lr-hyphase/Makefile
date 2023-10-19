IMAGE_NAME = lr-hifiphase
VERSION = 0.7.2

TAG1 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):$(VERSION)
#TAG2 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):latest

all: | build push

build:
    docker build -t $(TAG1) .

build_no_cache:
    docker build --no-cache -t $(TAG1) .

push:
    docker push $(TAG1)
