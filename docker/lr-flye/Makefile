VERSION = 2.8.3  # should match the Flye version number
TAG1 = quay.io/broad-long-read-pipelines/lr-flye:$(VERSION)
TAG2 = quay.io/broad-long-read-pipelines/lr-flye:latest
TAG3 = us.gcr.io/broad-dsp-lrma/lr-flye:$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/lr-flye:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
	docker push $(TAG3)
	docker push $(TAG4)
