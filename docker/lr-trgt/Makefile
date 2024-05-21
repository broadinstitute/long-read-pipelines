VERSION = 1.0.0
TAG1 = us.gcr.io/broad-dsp-lrma/lr-trgt:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-trgt:latest

all: build push

build:
        docker build -t $(TAG1) -t $(TAG2) .

build_no_cache:
        docker build --no-cache -t $(TAG1) -t $(TAG2) .

push:
        docker push $(TAG1)
        docker push $(TAG2)
