IMAGE_NAME = lr-c3poa
VERSION = 2.2.3 # This should match the version number from the C3POa repo

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
