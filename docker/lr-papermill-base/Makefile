VERSION = 2.5.0
TAG1 = us.gcr.io/broad-dsp-lrma/lr-papermill-base:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-papermill-base:latest

all: build

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
